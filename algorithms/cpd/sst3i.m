function [result]=sst3i(xi,M,K,l,p,q,bMu,nMu,varargin)
% Finds the points of change of x's dynamics incrementally. The output is
% returned as a NaN until there is enough data (q+M-1 points) to start the
% change discovery process.
% . Based on the following papers:
% 1. "Change Point Detection in Time Series by Means of
%    Singular Spectrum Analysis" V.G. Moskvina and A.A. Zhigljavsky
% 2. Application of the singular-spectrum analysis to change-point
%    detection in time series by the same authors (this is the main paper)
% 3. Software from http://www.cardiff.ac.uk/maths/subsites/stats/changepoint/
%    which is in the helpers folder of this package
%
% xi = NaN means reset. Otherwise it is the current point of the time
%     series.
% M = the width of the dynamic window (this is M in the paper)
% K = the number of windows to cosider in the past (this K in the paper)
%     here K>=M+1 because K=N-M+1 and N>=2M ==> K+M-1>=2M ==> K>=M+1
% p = The beginning of the test window (p in paper) default is K+1=N-M+2.
%     give -1 to keep the default
% q = the number of validation samples  (q in paper) default is K+2M-1=N+M.
%     give -1 to keep the default
% bMu = the beginning of the samples used to measure the mean distance
%       (default 0). will be used only if the optional parameter normalize
%       was nonzero which is the default but can be changed using optional
%       parameters
% nMu = the number of vectors in the mean calculatation (default K).will be used only if the optional parameter normalize
%       was nonzero which is the default but can be changed using optional
%       parameters
%
%
% For more information please consult the following publications: 
% ===============================================================
% V. Moskvina and A. Zhigljavsky. (2003) Change-point detection algorithm
% based on the singular-spectrum analysis, Communication in Statistics: 
% Simulation & Computation,  32, 319-352
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%

% 
persistent x;
x=[x,xi];
epsilon=1e-12;

N=M+K-1;
T=numel(x);
if(nargin <4)    
    l=-1;   
else
    if(l==0)
        l=-1;
    end;
end


if(nargin <5)
    p=N-M+2;
else
    if(p<0)
        p=N-M+2;
    end
end
if(nargin <6)
    q=N+M;
else
    if(q<0)
        q=N+M;
    end
end
Q=q-p;

if(nargin <7)
    bMu=0;
else
    if(bMu<0)
        bMu=0;
    end
end
if(nargin <8)
    nMu=K;
else
    if(nMu<0)
        nMu=K;
    end
end

maxNeeded=q+M-1;
if T>maxNeeded
    x=x(end-maxNeeded+1:end);
    T=N;
elseif T<maxNeeded
    result=NaN;
    return;
end

center=1;
useNoBalance=0;
normalizeS=1;
doth=0;
th=8/sqrt(M*(q-M-p));
lFun=@untilLargestGab;
lFunCalcOnce=1;
lFunParams=[];
lFunType=2;         % default apply to the power
nArgs=size(varargin,2);
if(nArgs>0)
    if(mod(nArgs,2)~=0)
        error('The optional arguments must be in the form name,value so they must be even!!!');
    end
    for i=1:2:nArgs
        switch(lower(varargin{i}))            
            case {'center'}                      % if nonzero then the output is cetnered which means that at time t the score is calculated using the past and the future not all in the future which is the default of this algorithm (default:center=1)
                center=varargin{i+1};            
            case 'nobalance'                     % if nonzero, the eigen vectors will be calculated without balance
                useNoBalance=varargin{i+1};            
            case {'n','normalize'}               % if nonzero, the statistic is normalized by division by the average distance in the test window
                normalizeS=varargin{i+1};            
            case {'ath','autothreshold'}         % if nonzero, the cumsum statistic (second output) will be thresholded automatically using a value from the software we are based on
                doth=varargin{i+1};            
            case {'th','threshold'}              % if the threshold to be used if authothreshold is nonzero. default is 4/sqrt(M*(q-M-p))
                th=varargin{i+1};                        
            case 'lfun'                         % useful only if l was passed <=0 and it gives the function used to calculate l. default is untilLargestGab. this function must have the same signature as untilLargestGab
                lFun=varargin{i+1};
            case 'lfunparams'                   % any parameters to be passed to lFun
                lFunParams=varargin{i+1};
            case 'lfuntype'                     % if 0 the lFun will be applied to the eigen values. otherwise they will be applied to eigen^lfuntype
                lFunType=varargin{i+1};
            case 'lfuncalconce'                  % if 0 then l will be recalculated every  iteration using the hankel matrix. if >0 then the l value will be calculated once using this number of points or the whole time series if it was -1                
                lFunCalcOnce=varargin{i+1};                        
            otherwise
                error('unknown command: %s',varargin{i});
        end
    end
end

% check inputs
if K<M+1
    error('K must be >= M+1!!!');
end
% if(p>=K)
% 	warning('Test samples should start after the training samples (p>=K)');
% end
if(q<N)
	error('Test samples must be at least exceeding training samples (q>=N)');
end
% if(Q>M)
% 	warning('Q(q-p) should be <= M');
% end
if(l>min(K,M))
	error('we cannot use eigen vectors more than the lowest of the hankel matrix dimensions');
end
if(M>N/2)
	error('Lag parameter (M) must be at most half the window width (N=M+K-1)');
end
if(bMu<0 || bMu>M)
	error('bMu cannot be negative or more than M');
end
if(nMu<0)
	error('nMu cannot be negative');
end

% for simplicity we do not try to calculate l automatically in the
% incremental version because there is no enough data in the beginning
% if(l<0)
%     if abs(lFunCalcOnce)<eps        
%     elseif lFunCalcOnce<0
%         X=zeros(M,T);
%         for j=1:T-M+1
%             X(:,j)=x(j:j+M-1);
%         end    
%         [U,S,V]=svd(X);
%         if abs(lFunType)<eps
%             l=lFun(diag(S),lFunParams);
%         else
%             l=lFun((diag(S).^2).^lFunType,lFunParams);
%         end
%     else
%         X=zeros(M,lFunCalcOnce);
%         for j=1:lFunCalcOnce-M+1
%             X(:,j)=x(j:j+M-1);
%         end       
%         [U,S,V]=svd(X);
%         if abs(lFunType)<eps
%             l=lFun(diag(S),lFunParams);
%         else
%             l=lFun((diag(S).^2).^lFunType,lFunParams);
%         end
%     end   
%     clear U;
%     clear V;
% end

% we do not need to do a shift in the incremental version
%outShift=1+(center~=0)*(N-1);

X=zeros(M,K);

%for n=0:T-N-1+K-q
n=0;
    % construct the Herkel Matrix of the future
    for j=1:K
        X(:,j)=x(n+j:n+j+M-1);
    end
    % calculate eigen vectors of R = XX'
    if useNoBalance
        [v,d]=eig((X*X'),'nobalance');
    else
        [v,d]=eig((X*X'));
    end
    % sort the eigen vectors
    d2=diag(d);
    [d2,di]=sort(d2,'descend');
    v=v(:,di);    
    % calculate l if needed
    if(l<0)
        if abs(lFunType)<eps
            l2=lFun(d2,lFunParams);
        else
            l2=lFun(abs(d2).^lFunType,lFunParams);
        end        
        U=v(:,1:l2);
    else
        % keep the top l eigen vectors
        U=v(:,1:l);        
    end
    
    U2=(U*U');
    % calculate the average distance between vectors of the test matrix and
    % the subspace spanned by U
    D=0;    
    for j=p+1:q
        Xj=x(n+j:n+j+M-1)';
        D=D+ Xj'*Xj-Xj'*U2*Xj;
    end
    D=D./(M*(q-p));    
    
    % calculate normalization constant if needed
    if(abs(normalizeS)<eps)
        mu=1;
    else
        mu=0;
        for j=bMu+1:bMu+nMu
            Xj=x(n+j:n+j+M-1)';
            mu=mu+ Xj'*Xj-Xj'*U2*Xj;
        end
        mu=mu./(M*(nMu));        
    end    
    if(abs(mu)<epsilon && abs(D)<epsilon)
        result=0;    
    elseif abs(mu)<epsilon
        result=0; %result(n+outShift-1);
    else        
        result=D./mu;    
    end    
%end

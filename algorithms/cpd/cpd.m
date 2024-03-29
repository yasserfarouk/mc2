% This file is a part of the MC2 toolbox developed by Y. Mohammand and T. Nishida.
%Please do not remove this comment
%
% Using this file is governed by the license of MC2 which you can find in LICENSE.md
% 
% You can find more information about this toolbox here:
% - Yasser Mohammad and Toyoaki Nishida, "MC2: An Integrated Toolbox for Change, Causality, 
%   and Motif Discovery", 29th International Conference on Industrial, Engineering & 
%   Other Applications of Applied Intelligent Systems (IEA/AIE) 2016, pp. 128 -- 141.
% - Yasser Mohammad and Toyoaki Nishida, "Data Mining for Social Robotics", Springer 2016.
%

function [result,chgLocs,cumsum,chgLocsCumsum]=cpd(x,M,K,l,varargin)
% Runs several variants of SSA based change point discovery. The exact
% variant can be selected by setting optional parameters
%
% Inputs:
% ======
% x = the input vector 
% M = the width of the dynamic window (this is M in the paper)
% K = the number of windows to cosider in the past (this K in the paper)
%     here K>=M+1 because K=N-M+1 and N>=2M ==> K+M-1>=2M ==> K>=M+1
% varargin = see the switch statement in the beginning of the routine
%
% Outputs:
% =======
% result    A time series of the same length as the input x but with change
%           score at every point
% chgLocs   [Optional] The localized change locations of the input
% cumsum    [Optional] A cumsum like statistic based on result
% chgLocsCumsum [Optional] localized change locaitons in the cumsum
%                           statistic output
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, On Comparing SSA-based Change Point
% Discovery Algorithms, IEEE SII 2011, 938-945 
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%

% 
% 
global savedUs;
global d2s;
global reuseUs;
global savedTo;

epsilon=1e-30;

if(size(x,2)==1)
    x=x';
end

maxOutput=inf;
minOutput=0;
symmetricNormalization=[];
T=numel(x);
N=M+K-1;
w=2*N+1;
nKf=K;
nKp=K;
gff=0;
gp=0;
if nargin>3 && ischar(l) && (mod(length(varargin),2)~=0)
    varargin=[l,varargin];    
    l=-1;
end
%nK=K;
if(nargin <4)    
    l=-1;   
else
    if(l==0)
        l=-1;
    end;
end

center=1;
useNoBalance=0;
normalizeS=1;
scaleOutput=1;
useSVD=1;
calcHEvery=1;
hType='first';
doThCS=0;
th=[];
reuseUs=1;

lMax=min(M,K);
lMin=1;
lFun=@untilLargestGab;
nPoints4H=0;
lFunParams=[];
lFunType=0;         % default apply to the power

doPost=0;
pFun=@rsstpost;
pFunParams=1;

locFunCS=@findLocsCS;
locFunCSParams=[];

locFun=@findLocsTh;
locFunParams=[]; % the default is to calculate the mean point of change

distFun=@distAvgVHDist;
distFunParams=[];

direction ='forward';
showProgress=0;
nArgs=size(varargin,2);
if(nArgs==1 && iscell(varargin{1}))
    varargin=varargin{1};
    nArgs=size(varargin,2);
end
if(nArgs>0)
    if(mod(nArgs,2)~=0)
        error('The optional arguments must be in the form name,value so they must be even!!!');
    end
    for i=1:2:nArgs
        switch(lower(varargin{i}))
            case {'sst3defaults'}
                if(varargin{i+1})
                    nKf=2*M-2;
                    gff=2-M;
                    lFunType=2;
                    useSVD=0;
                    nPoints4H=0;
                    calcHEvery=1;
                    normalizeS=1;
                    distFun=@distAvgVHDist;
                    scaleOutput=0;
                end
            case {'sstdefaults'}        %may require correction
                if(varargin{i+1})                    
                    lFunType=0;
                    useSVD=0;
                    nPoints4H=-1;
                    calcHEvery=1;
                    normalizeS=0;
                    distFun=@distAvgVHDist;
                    scaleOutput=0;
                end
            case {'rsstdefaults'}
                if(varargin{i+1})                                        
                    nPoints4H=0;
                    calcHEvery=1;
                    normalizeS=1;
                    distFun=@distWeightedEigDist;
                    scaleOutput=0;
                    doPost=1;
                end
            case {'symnorm','symmetricnormalization'} % only of value if normalization is active. If nonzero then normalization will be done by (future-past)/(future+past). By default it is nonzero of calcHEvery<=0
                symmetricNormalization=(varargin{i+1});   
            case {'maxputout'}              % maximum output score (scores above that are truncated to it) default inf
                maxOutput=(varargin{i+1});   
            case {'minputout'}              % minimum output score (scores below that are increased to it) default 0
                minOutput=(varargin{i+1});   
            case {'scaleoutput'}
                scaleOutput=(varargin{i+1});   
            case {'dir','direction'}        % direction of applying the algorithm which is either 'forward' or 'backward'. default is 'forward'
                direction=(varargin{i+1});   
            case {'reuseus','fast'}         % if nonzero the hankel matrix will be reused whenever possible
                reuseUs=(varargin{i+1});   
            case {'lmax','maxl'}            % the maximum number of eigen vectors allowed in any hankel matrix. a negative value means ignore this limit (default -1)
                lMax=(varargin{i+1});   
            case {'lmin','minl'}            % the maximum number of eigen vectors allowed in any hankel matrix. a negative value means ignore this limit (default -1)
                lMin=(varargin{i+1});   
            case {'kf','nkf'}              % the number of vectors to be used from the future. Default is K
                nKf=(varargin{i+1});    
            case {'kp','nkp'}              % the number of vectors to be used from the past. Default is K
                nKp=(varargin{i+1});    
            case {'gf','futureskip'}             % the number of points to skip in the future when searching for a change (default 0)
                gff=(varargin{i+1});    
            case {'gp','pastskip'}               % the number of points to skip in the past when searching for a change (default 0)
                gp=(varargin{i+1});    
            case 'useeig'                        % if nonzero then eigen values will be used for alalyzing hankel matrix otherwise SVD
                useSVD=(varargin{i+1});    
            case 'usesvd'                        % if nonzero then SVD values will be used for alalyzing hankel matrix otherwise eigen values
                useSVD=varargin{i+1};    
            case 'dopost'                        % if nonzero then post processing is applied
                doPost=varargin{i+1};    
            case {'dfun','distfun'}             % the distance function to be used for comparing the past and the future. See distAvgHVDist for the declaration
                distFun=varargin{i+1};
            case {'dfunname','distfunname'}             % the distance function to be used for comparing the past and the future. See distAvgHVDist for the declaration
                distFun=str2func(varargin{i+1});
            case {'dfparams','dparams','dfunparams'} % any parameters to be passed to pFun
                distFunParams=varargin{i+1};            
            case {'pfun','postfun'}             % useful only if doPost was passed as nonzero. this function receives the result and M and K and may receive other parameters if pFunParams was given
                pFun=varargin{i+1};
            case {'pfparams','pparams','pfunparams'} % any parameters to be passed to pFun
                pFunParams=varargin{i+1};            
            case {'locfunth','locfun'}                 % useful only if the second output is required. this function receives the result and any parameters passed to it using locFunParams and returns a vector of locations of change
                locFun=varargin{i+1};
            case {'locfthparams','locthparams','locfunthparams','locfparams','locparams','locfunparams'} % any parameters to be passed to locfun. if th was passed then it will be passed to locFun before these parameters
                locFunParams=varargin{i+1};            
            case {'th','threshold','pchange'}       % the threshold to be used for localizing changes in the time series. will be passed as first locFunCSParams if specified
                th=varargin{i+1};            
            case {'locfuncs'}             % useful only if the fourth output is required. this function receives the cumsum and any parameters passed to it using locFunParams and returns a vector of locations of change
                locFunCS=varargin{i+1};
            case {'locfcsparams','loccsparams','locfuncsparams'} % any parameters to be passed to locfuncs
                locFunCSParams=varargin{i+1};            
            
            case {'center'}                      % if nonzero then the output is cetnered which means that at time t the score is calculated using the past and the future not all in the future which is the default of this algorithm (default:center=1)
                center=varargin{i+1};            
            case 'nobalance'                     % if nonzero, the eigen vectors will be calculated without balance
                useNoBalance=varargin{i+1};            
            case {'normalize'}               % if nonzero, the statistic is normalized by division by the average distance in the test window
                normalizeS=varargin{i+1};            
            case {'ath','autothreshold'}         % if nonzero, the cumsum statistic (second output) will be thresholded automatically using a value from the software we are based on
                doThCS=varargin{i+1};            
            case {'csth','csthreshold'}         % if the threshold to be used if authothreshold is nonzero. default is 4/sqrt(M*(q-M-p))
                csth=varargin{i+1}; 
            case {'calchankelevery','lfuncalconce'}  % if ==0 then the hyperplane used to compare all points is calculated once in the beginning using the whole time series. if < 0 it will be calculated once using this number of subsequences selected from the time series based on the randomH value. if >0 then it will be calculated every specified number of points (default 1)
                calcHEvery=(varargin{i+1});    
            case 'htype'                       % if calchankelevery<0. if 'first' then first subsequences will be used. if 'rand' then random subsequences from the series will be used. if 'last' then last subsequences will be used. if 'equidist' then the vectors will be selected to cover the whole series with equal time spans between them. if 'center' the center of the series will be used
                hType=(varargin{i+1});      % if calchankelevery>0. it affects the same way the evaluation of the hankel matrix     
            %case 'cmptype'                     % the way to select the vectors to compare past and future. 'first', 'last','rand','equidist' same meaning as hType,'around' means the past is taken last and future first, 'far' means that the past is taken 'first' and the future 'last'
             %   cType=(varargin{i+1}==0);      
            case 'lfun'                         % useful only if l was passed <=0 and it gives the function used to calculate l. default is untilLargestGab. this function must have the same signature as untilLargestGab
                lFun=varargin{i+1};
            case 'lfunparams'                   % any parameters to be passed to lFun
                lFunParams=varargin{i+1};
            case 'lfuntype'                     % if 0 the lFun will be applied to the eigen values. otherwise they will be applied to eigen^lfuntype
                lFunType=varargin{i+1};
            case {'npts4h','npointsforhankel'}  % if 0 then l will be recalculated every  iteration using the hankel matrix. if >0 then the l value will be calculated once using this number of points or the whole time series if it was -1
                nPoints4H=varargin{i+1};       
            case 'showprogress'
                showProgress=varargin{i+1};       
            otherwise
                error('unknown command: %s',varargin{i});
        end
    end
end

if isempty(symmetricNormalization)
    symmetricNormalization=abs(calcHEvery)<eps;
end

wMax=max(w,N);              % this is the end of the window containing both the henkel matrix and the future and past vectors
mdPnt=int32(floor(w/2));
outShift=1+int32((center~=0))*(mdPnt-1);

p=mdPnt+gff;
q=p+nKf;
bMu=mdPnt-gp-nKp-M+1;
nMu=nKp;

csth=8/sqrt(double(M*(q-M-p)));     % todo find a better threshold

% check inputs
if K<M+1
    warning('K should be >= M+1!!!');
end
% if(p>=K)
% 	warning('Test samples should start after the training samples (p>=K)');
% end
% if(q<N)
% 	error('Test samples must be at least exceeding training samples (q>=N)');
% end
% if(Q>M)
% 	warning('Q(q-p) should be <= M');
% end
if(l>min(K,M))
	error('we cannot use eigen vectors more than the lowest of the hankel matrix dimensions');
end
if(M>N/2)
	error('K must be > M');
end
% if(bMu<0 || bMu>M)
% 	error('bMu cannot be negative or more than M');
% end
% if(nMu<0)
% 	error('nMu cannot be negative');
% end


if showProgress
    fprintf('starting CPD on %d points\n',T);
end

U=[];
U2=[];
if(l<0)
    if abs(nPoints4H)<eps        
    else
        if nPoints4H<0          % if this is changed then check calcHEvery block beneath it.
            X=zeros(M,T);
            for j=1:T-M+1
                X(:,j)=x(j:j+M-1);
            end                
        else
            X=zeros(M,nPoints4H);
            for j=1:nPoints4H-M+1
                X(:,j)=x(j:j+M-1);
            end                   
        end
        if( useSVD~=0)
            [U,S,V]=svd(X);
            d2=diag(S).^2;
        else
            if useNoBalance
                [v,d]=eig((X*X'),'nobalance');
            else
                [v,d]=eig((X*X'));
            end
            % sort the eigen vectors
            d2=diag(d);
            [d2,di]=sort(d2,'descend');
            U=v(:,di);    
        end
        if abs(lFunType)<eps
            l=lFun(d2,lFunParams);
        else
            l=lFun((d2).^lFunType,lFunParams);
        end        
        if(lMax>0)
            l=min(l,lMax);
        end
        if(lMin>0)
            l=max(l,lMin);
        end
        clear V;
    end            
end

if(abs(calcHEvery)<eps)
    if(nPoints4H<0 && ~isempty(U)) % U is already calculated then do nothing        
        U=U(:,l);
    else
        X=zeros(M,T);
        for j=1:T-M+1
            X(:,j)=x(j:j+M-1);
        end    
        if( useSVD~=0)
            [v,d,dummy]=svd(X);        
            d2=diag(d).^2;            
        else
            if useNoBalance
                [v,d]=eig((X*X'),'nobalance');
            else
                [v,d]=eig((X*X'));
            end
            % sort the eigen vectors
            d2=diag(d);
            [d2,di]=sort(d2,'descend');
            v=v(:,di);    
        end
        % calculate l if needed
        if(l<0)
            if abs(lFunType)<eps
                l2=lFun(d2,lFunParams);
            else
                l2=lFun(abs(d2).^lFunType,lFunParams);
            end        
            if(lMax>0)
                l2=min(l2,lMax);
            end
            if(lMin>0)
                l2=max(l2,lMin);
            end
            U=v(:,1:l2);
        else
            % keep the top l eigen vectors
            U=v(:,1:l);        
        end            
    end
    U2=U*U';    
elseif calcHEvery<0
    nVects=-calcHEvery;
    nPts=nVects+M-1;
    X=zeros(M,nPts);
    switch(hType)
        case 'first'            
            for j=1:nVects
                X(:,j)=x(j:j+M-1);
            end                
        case 'last'
            for j=T-nPts+1:T
                X(:,j)=x(j:j+M-1);
            end
        case 'rand'
            rndV=(T-M).*rand(nVects,1)+1;
            for j=1:nVects                
                X(:,j)=x(rndV(j):rndV(j)+M-1);
            end                
        case 'equidist'
            rndV=1:(T-M+1)/(nVects-1):T-M+1;
            for j=1:nVects                
                X(:,j)=x(rndV(j):rndV(j)+M-1);
            end
        case 'center'
            beg=int32(floor(T/2-nVects/2));
            for j=beg:beg+nVects-1
                X(:,j)=x(j:j+M-1);
            end                
        otherwise
            error(sprintf('the hankel matrix type %s is not supported',hType));
    end    
    
    if( useSVD~=0)
        [v,d,dummy]=svd(X);        
        d2=diag(d).^2;        
    else
        if useNoBalance
            [v,d]=eig((X*X'),'nobalance');
        else
            [v,d]=eig((X*X'));
        end
        % sort the eigen vectors
        d2=diag(d);
        [d2,di]=sort(d2,'descend');
        v=v(:,di);    
    end
    % calculate l if needed
    if(l<0)
        if abs(lFunType)<eps
            l2=lFun(d2,lFunParams);
        else
            l2=lFun(abs(d2).^lFunType,lFunParams);
        end        
        if(lMax>0)
            l2=min(l2,lMax);
        end
        if(lMin>0)
            l2=max(l2,lMin);
        end
        U=v(:,1:l2);
    else
        % keep the top l eigen vectors
        U=v(:,1:l);        
    end    
    U2=(U*U');
    wMax=w;         % we will calculate no matrices in the loop
else
    % do nothing it will be calculated in the loop
    switch(hType)
        case 'first'           
            assert(N<=w,sprintf('you cannot use hType=%s while N(%d) is larger than w(%d)',hType,N,w));
            vecIndices=1:K;
        case 'last'
            assert(N<=w,sprintf('you cannot use hType=%s while N(%d) is larger than w(%d)',hType,N,w));
            vecIndices=w-M-K+2:w-M+1;
        case 'rand'
            vecIndices=int32((w-M).*rand(1,K)+1);
        case 'equidist'
            vecIndices=1:int32((w-M)/(K-1)):w-M+1;
        case 'center'
            b=int32(floor(N/2-K/2));
            vecIndices=b:b+K-1;
        otherwise
            error(sprintf('the hankel matrix type %s is not supported',hType));
    end
end

savedUs=cell(T,1);      % may be used to save hankel matrix eigen vectors if possible
d2s=cell(T,1);
savedTo=zeros(T,1);
% if reuseUs
%     savedUs=cell(T,1);      % may be used to save hankel matrix eigen vectors if possible
%     d2s=cell(T,1);
%     savedTo=cell(T,1);
% else
%     savedUs=[];
%     d2s=[];
%     savedTo=[];
% end

if(strcmpi(direction,'backward'))
    x=x(T:-1:1);
end

X=zeros(M,K);
nxtHCalc=0;
result=zeros(1,T);
sstart=tic;
showProgressEvery=floor(((T-wMax))*(showProgress/100));
if showProgress
ssss=sprintf('%03.2f,remaining time:%03.2f minutes',999,99);
fprintf(ssss);
end

for t=0:T-wMax%     T-N-1+K-q 
    if showProgress
        fractionCompleted=t./(T-wMax);
        if ~mod(t,showProgressEvery) && t>0
            tm=toc(sstart);
            fprintf(repmat('\b',1,length(ssss)));
            ssss=sprintf('%03.2f,remaining time:%03.2fminutes',fractionCompleted*100,(tm/60)*((T-wMax-t)/t));            
            fprintf(ssss);
        end
    end
    if(abs(calcHEvery)<eps)
    elseif calcHEvery>0
        if(t>=nxtHCalc)
            nxtHCalc=nxtHCalc+calcHEvery;
            if(~reuseUs ||( isempty(savedUs{t+vecIndices(1)}) || savedTo(t+vecIndices(1))~=K))
                % construct the Herkel Matrix of the past
                for j=1:K
                    X(:,j)=x(t+vecIndices(j):t+vecIndices(j)+M-1);
                end
                % calculate eigen vectors of R = XX'
                if( useSVD~=0)
                    [v,d,dummy]=svd(X);        
                    d2=diag(d);
                    d2=d2.*d2;
                else
                    if useNoBalance
                        [v,d]=eig((X*X'),'nobalance');
                    else
                        [v,d]=eig((X*X'));
                    end
                    % sort the eigen vectors
                    d2=diag(d);
                    [d2,di]=sort(d2,'descend');
                    v=v(:,di);    
                end
                % calculate l if needed
                if(l<0)
                    if abs(lFunType)<eps
                        l2=lFun(d2,lFunParams);
                    else
                        l2=lFun(abs(d2).^double(lFunType),lFunParams);
                    end        
                    if(lMax>0)
                        l2=min(l2,lMax);
                    end
                    if(lMin>0)
                        l2=max(l2,lMin);
                    end
                    U=v(:,1:l2);
                else
                    % keep the top l eigen vectors
                    U=v(:,1:l);        
                end    
                savedUs{t+vecIndices(1)}=U;
                savedTo=K;
                d2s{t+1}=d2;
            else
                U=savedUs{t+vecIndices(1)};
   %             d2=d2s{t+1};
            end
            U2=(U*U');
            
        end
    end
    % calculate the average distance between vectors of the test matrix and
    % the subspace spanned by U
    D=distFun(x,t,M,K,U, p,q,l,lMax,lMin,lFun, lFunType, lFunParams,useSVD, useNoBalance,U2,distFunParams);
    
    % calculate normalization constant if needed
    if(abs(normalizeS)<eps)        
        result(t+outShift)=D;
    else
        mu=distFun(x,t,M,K,U,bMu,bMu+nMu,l,lMax,lMin,lFun, lFunType, lFunParams,useSVD, useNoBalance,U2,distFunParams);        
        if(symmetricNormalization)
            if(abs(mu)<epsilon && abs(D)<epsilon)
                result(t+outShift) =0;    
            else
                mx=max([mu,D]); mn=min([mu,D]);
                if abs(mn)<epsilon
                    result(t+outShift)=0; %min(D,max(result(1:t+outShift-1)));
                else        
                    result(t+outShift)=(mx-mn)./(mx+mn);    
                end
            end
        else
            if(abs(mu)<epsilon && abs(D)<epsilon)
                result(t+outShift) =0;    
            elseif abs(mu)<epsilon
                result(t+outShift)=0; %min(D,max(result(1:t+outShift-1)));
            else        
                result(t+outShift)=D./mu;    
            end
        end
    end    
    
    %assert(~isnan(result(t+outShift)) && ~isinf(result(t+outShift)));
end

if(strcmpi(direction,'backward'))
    result=result(T:-1:1);
end 
result(result>maxOutput)=maxOutput;
result(result<minOutput)=minOutput;
if scaleOutput
    if(max(result)>epsilon)
        result=result./max(result);
    end
    result(result<epsilon)=0;
end
%result(result>1)=1;

%do post processing using same M and K
if doPost~=0
    result=pFun(result,M,K,pFunParams);
end

%extra thining operation to keep only local minima
if nargout>1
    if(~isempty(th))        
        chgLocs = locFun(result,th,deal(locFunParams));
    else
        chgLocs = locFun(result,deal(locFunParams));
    end
end

% calculate cumsum if needed
if(nargout>2)
    cumsum=result;
    %a=mean(result);
    for i=1:T-1
        %cumsum(i+1)=cumsum(i)+result(i+1)-result(i)-1./sqrt(M*(q-p));
        cumsum(i+1)=max(0,cumsum(i)+result(i+1)-1);
    end
    for i=1:T
        if(doThCS)
            if(cumsum(i)<csth)
                cumsum(i)=0;
            else
                cumsum(i)=1;
            end
        end
    end
end

if nargout>3
    if(doThCS)
        chgLocsCumsum = locFunTh(cumsum,locFunCSParams);
    else        
        chgLocsCumsum = locFunCS(cumsum,locFunThParams);
    end
end
if showProgress
    fprintf('\n');
end
clear savedUs;
clear d2s;
clear reuseUs;
clear savedTo;
end

function locs=findLocs(r,varargin)
    r=thin(r);
    locs=find(abs(r)>0.000001);
end

function d=distAvgVHDistNorm(x,t,M,K,U,p,q,l,lMax,lMin,lFun, lFunType, lFunParams,useSVD, useNoBalance,U2,varargin)
%varargin [optional] U2    
    d=0.0;    
    for j=p+1:q
        Xj=x(t+j:t+j+M-1)'; 
        Xj=Xj./norm(Xj);
        d=d+ abs(1-Xj'*U2*Xj);
    end
    d=d./double((q-p));    
end

function d=distAvgVHDist(x,t,M,K,U,p,q,l,lMax,lMin,lFun, lFunType, lFunParams,useSVD, useNoBalance,U2,varargin)
%varargin [optional] U2    
    d=0.0;    
    for j=p+1:q
        Xj=x(t+j:t+j+M-1)';            
        d=d+ abs(Xj'*Xj-Xj'*U2*Xj);
    end
    d=d./double(M*(q-p));    
end

function d=distAngBetweenSubspaces(x,t,M,K,U,p,q,l,lMax,lMin,lFun, lFunType, lFunParams,useSVD, useNoBalance,U2,varargin)
global savedUs;
global d2s;
global savedTo;
global reuseUs;
if(~reuseUs ||( isempty(savedUs{t+p+1}) || (savedTo(t+p+1)~=(q-p)) ))
    X=zeros(M,q-p);
    for j=p+1:q
        X(:,j)=x(t+j:t+j+M-1)';                    
    end

    if( useSVD~=0)
        [v,d,dummy]=svd(X);        
        d2=diag(d).^2;        
    else
        if useNoBalance
            [v,d]=eig((X*X'),'nobalance');
        else
            [v,d]=eig((X*X'));
        end
        % sort the eigen vectors
        d2=diag(d);
        [d2,di]=sort(d2,'descend');
        v=v(:,di);    
    end
    % calculate l if needed
    if(l<0)
        if abs(lFunType)<eps
            l2=lFun(d2,lFunParams);
        else
            l2=lFun(abs(d2).^lFunType,lFunParams);
        end        
        if(lMax>0)
            l2=min(l2,lMax);
        end
        if(lMin>0)
            l2=max(l2,lMin);
        end
        Unew=v(:,1:l2);
    else
        % keep the top l eigen vectors
        Unew=v(:,1:l);        
    end            
    savedUs{t+p+1}=Unew;
    d2s{t+p+1}=d2;        
    savedTo(t+p+1)=q-p;
else
    Unew=savedUs{t+p+1};    
end
    d=1-abs(cos(subspace(U,Unew)));    
end

function d=distAvgEigDist(x,t,M,K,U,p,q,l,lMax,lMin,lFun, lFunType, lFunParams,useSVD, useNoBalance,U2,varargin)
global savedUs;
global d2s;
global savedTo;
global reuseUs;
if(~reuseUs ||( isempty(savedUs{t+p+1}) || (savedTo(t+p+1)~=(q-p)) ))
    X=zeros(M,q-p);
    for j=p+1:q
        X(:,j)=x(t+j:t+j+M-1)';                    
    end

    if( useSVD~=0)
        [v,d,dummy]=svd(X);        
        d2=diag(d).^2;        
    else
        if useNoBalance
            [v,d]=eig((X*X'),'nobalance');
        else
            [v,d]=eig((X*X'));
        end
        % sort the eigen vectors
        d2=diag(d);
        [d2,di]=sort(d2,'descend');
        v=v(:,di);    
    end
    % calculate l if needed
    if(l<0)
        if abs(lFunType)<eps
            l2=lFun(d2,lFunParams);
        else
            l2=lFun(abs(d2).^lFunType,lFunParams);
        end        
        if(lMax>0)
            l2=min(l2,lMax);
        end
        if(lMin>0)
            l2=max(l2,lMin);
        end
        Unew=v(:,1:l2);
    else
        % keep the top l eigen vectors
        Unew=v(:,1:l);        
    end            
    savedUs{t+p+1}=Unew;
    d2s{t+p+1}=d2;        
    savedTo(t+p+1)=q-p;
else
    Unew=savedUs{t+p+1};
end
    d=0;
    l=size(Unew,2);
    for i=1:l
        d=d+ abs(Unew(:,i)'*Unew(:,i)-Unew(:,i)'*U2*Unew(:,i));
    end
    d=d/double(l);
end

function d=distWeightedEigDist(x,t,M,K,U,p,q,l,lMax,lMin,lFun, lFunType, lFunParams,useSVD, useNoBalance,U2,varargin)
global savedUs
global d2s
global reuseUs
global savedTo
    if(~reuseUs ||( isempty(savedUs{t+p+1}) || (savedTo(t+p+1)~=(q-p)) ))
        X=zeros(M,q-p);
        for j=p+1:q
            X(:,j)=x(t+j:t+j+M-1)';                    
        end

        if( useSVD~=0)
            [v,d,dummy]=svd(X);        
            d2=diag(d).^2;        
        else
            if useNoBalance
                [v,d]=eig((X*X'),'nobalance');
            else
                [v,d]=eig((X*X'));
            end
            % sort the eigen vectors
            d2=diag(d);
            [d2,di]=sort(d2,'descend');
            v=v(:,di);    
        end
        % calculate l if needed
        if(l<0)
            if abs(lFunType)<eps
                l2=lFun(d2,lFunParams);
            else
                l2=lFun(abs(d2).^lFunType,lFunParams);
            end        
            if(lMax>0)
                l2=min(l2,lMax);
            end
            if(lMin>0)
                l2=max(l2,lMin);
            end
            Unew=v(:,1:l2);
        else
            % keep the top l eigen vectors
            Unew=v(:,1:l);        
        end        
        savedUs{t+p+1}=Unew;
        d2s{t+p+1}=d2;        
        savedTo(t+p+1)=q-p;
    else
        Unew=savedUs{t+p+1};
        d2=d2s{t+p+1};
    end
    d=0;
    l=size(Unew,2);    
    s=0;
    d2=d2./max(d2);
    for i=1:l
        d=d+ (d2(i).^2)*(Unew(:,i)'*Unew(:,i)-Unew(:,i)'*U2*Unew(:,i));        
        s=s+d2(i).^2;
    end
    if(abs(s)<eps)
        d=0;
    else
        d=d/double(s);
    end
end

function locs=findLocsCS(r,varargin)
    r=thin(r);
    locs=find(abs(r)>0.000001);
end


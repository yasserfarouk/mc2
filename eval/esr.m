function q=esr(x,t,nPoints,scalingFun,param)
% finds the Equivalence Sampling Rate between x and t pdfs given that t
% represents true probability of occurrence of some event.
%
% ESR is the probability that if I sample from x the sample will be within
% nPoints from some nonzero sample in t
%
% x         output of change point algorithm. If x is n*m matrix then it is
%           treated as m outputs of m algorithms with a time series of
%           length n and nPoints must be n points. This means the each
%           algorithm's output must be in a column not a row
% t        true change points (e.g. 1 at change points and 0 otherwise)
%           the value at each change point may be treated as weight and the
%           more this value is the more important this change point will be
%           in the result. This must be a vector
% nPoints   width around change points that is assumed OK. 
%
% scalingFun [Optional] a function name or handle that can be used to weight distances from true
%            value. it should take a vector of integers from -nPoints:+nPoints and
%            return a value from 0 to 1 that scores this delay. It may be
%            sign dependent but usually is not. The function can also
%            receive any arguments through param
%
%            The default scalingFun is a uniform distribution
%            You can use 'normalScaling' for a normal scaling function with
%            stddev equal to nPoints
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



    % set constants
    epsilon=1e-12;
    
    if nargin<4
        scalingFun=[];
    else
        if ischar(scalingFun)
            scalingFun=str2func(scalingFun);
        end
    end
    if nargin<5
        param=[];
    end
    % confirm input
    if ~isvector(t)
        error('t must be a vector');
    end
    if(size(x,1)==1) && (size(x,2)~=1)
        x=x';
    end


    % confirm and normalize
    n=size(x,1);
    m=size(x,2);
    if n ~= numel(t)
        error('The number of elements in t must be equal to the number of columns in x');
    end
    
    st=sum(t);    
    if(abs(st)<epsilon)
        error('the true pdf must contain at least a single nonzero value');
    end
    t=t./st;
    % compute
    q=zeros(m,1);
    qmax=esrhelper(t,t,n,nPoints,scalingFun,param);
    for i=1:m
        sx=sum(x(:,i)); 
        if(abs(sx)<epsilon)
            warning('ESR:AllZero','\t\tx(%d) is all zeros, we will report a zero for this signal',i);
            q(i)=0;
        else
            x(:,i)=x(:,i)./sx;
            q(i)=(0.5*esrhelper(x(:,i),t,n,nPoints,scalingFun,param) ...
                +0.5*esrhelper(t,x(:,i),n,nPoints,scalingFun,param))/qmax; 
        end        
    end        
end

function q=esrhelper(x,t,n,nPoints,scalingFun,param)    
    q=0;
    if(isempty(scalingFun))
        for i=1:n
            range=[max(1,i-nPoints) min(n,i+nPoints)];
            q=q+sum(t(range(1):range(2))).*x(i);
%            p=0;
%             for j=range(1):range(2)
%                 p=p+t(j);
%             end
%             p=sum(t(range(1):range(2)));
%             q=q+p.*x(i);
        end
    else        
        t=t(:);
        for i=1:n
            range=[max(1,i-nPoints) min(n,i+nPoints)];
            p=sum(t(range(1):range(2))'.*scalingFun([range(1)-i:range(2)-i],param));           
            p=p.*(range(2)-range(1)+1)./(sum(scalingFun([range(1)-i:range(2)-i])));
            q=q+p.*x(i);
        end 
    end
end

% example of a scaling function (normal with a standard deviation of 5
% points. it receives one onptional parameter (the variance). The mean is
% always zero
function w=normalScaling(d,param)
    if(nargin<2)
        stddev=5;
    else        
        stddev=param;
    end
    w=stddev.*normpdf(d);
end

function w=expScaling(d,param)
    if(nargin<2 || isempty(param))
        param=1;    
    end
    w=(exp(-abs(param.*d)));
end
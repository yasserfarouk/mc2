function [x,chgLocs,chgs,xNoiseless]=produceSingle(T,varargin)
% prodcues example 1-D time series of length T with random change locations.
%
% T             The length of the time series
% varargin      See the switch statement at the beginning of the code
%
%
% Outputs:
% ========
% x             The matrix  representing the time series (sum(nDim)*T)
% chgLocs       The locations of the changes
% chgs          An matrix of same size as x but with one at change
%               locations and zero everywhere else
% xNoiseless    x but without any noise or outliers added and without smoothing.
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

    pAR=0.1429;            % probability that the geneation is from an AR Model
    arParamRange=[-.5,.5];
    arModelRange=[2,20];
    
    minChangeDistance=ceil(T/100);    
    pchange=0.01;    
    range=[-inf, inf];
    scale2range=0;
    adjustEPA=1;    
    nArgs=size(varargin,2);
    noiseSigma=0;
    noiseLevel=0;
    pOutlier=0;
    outlierRange=[];
    smooth=0;
    if(nArgs>0)
        if(mod(nArgs,2)~=0)
            error('The optional arguments must be in the form name,value so they must be even!!!');
        end
        for i=1:2:nArgs
            switch(lower(varargin{i}))       
                case 'poutlier'                 % the probability of having outliers in the signal. default is 0
                    pOutlier=(varargin{i+1});
                case 'outlierrange'             % the range of values for outliers (default is range if scale2range is given and [-100,100] otherwise)
                    outlierRange=(varargin{i+1});
                case {'noiselevel','level'}     % The maximum of a uniform noise from -noiseLevel to +noiseLevel added to the signal. default is 0                    
                    noiseLevel=(varargin{i+1});
                case {'noisesigma','sigma'}     % the standard deviation or normal zero mean random noise to be added to the signal
                    noiseSigma=(varargin{i+1});
                case {'arparramrange','arprange'} % in order the minimum and maximum possible value of the parameter range. default is [-.5,.5] which ensures stability and not very large fluctuation
                    arParamRange=(varargin{i+1});
                case {'armodelrange','arrange'} % in order the minimum and maximum order of the AR model (default 2,minChangeDistance/10);
                    arModelRange=(varargin{i+1});
                case {'pc','pchange'}        % probability of a change. default is 0.01
                    pchange=(varargin{i+1});
                %case {'pchangenoparent','pcnp'}        % probability of a change with absolutely no reason. default 0
                 %   pchangeNoReason=(varargin{i+1});               
                case {'minChangeDistance','minchangedistance'}              % the minimum distance between any two allowable changes. default T/100. Notice that sometimes changes may happen faster than that if a time series had multiple parents (causes)
                    minChangeDistance=(varargin{i+1});                
                case {'pautoregressive','par'}         % the probability of using an autoregressive model
                    pAR=(varargin{i+1});     
                case {'range'}                         % the range to allow the output to be within. Default -Inf:Inf
                    range=(varargin{i+1});     
                case {'scale2range'}                   % if nonzero then the output will be scaled to take the whole range form range(1) to range (2). default: 0
                    scale2range=(varargin{i+1});     
                case {'adjusteverypatternalone'}       % if nonzero each pattern will be adjusted to the range alone otherwise, the whole signal will be adjusted in the end. default:1
                    adjustEPA=(varargin{i+1});     
                case {'smooth'}                        % if nonzero the output will be smoothed with minChangeDistance kernel. the default is no smoothing
                    smooth=(varargin{i+1});
                otherwise
                    error('Unknown argument: %s',varargin{i});
            end
        end
    end
    if(~adjustEPA)
        ss2range=scale2range;
        srange=range;
        range=[-inf, inf];
        scale2range=0;
    end
    chgs=zeros(1,T);
    chgLocs=[];
    % fill all time series that depends on nothing    
    
    % implement coupling    
    pchg=rand(1,T);
    l=1;
    while l<T-minChangeDistance
        if(pchg(l)<=pchange)
            chgLocs=[chgLocs;l];
            l=l+minChangeDistance;
        end
        l=l+1;
    end
    
    % do the changes in the correct time steps in all dimensions    
    % fill with fixed pattern
    if(rand(1,1)<pAR)
        param=(arModelRange(2)-arModelRange(1))*randi(1,1)+arModelRange(1);
        x=generateRandAR(T,param,range,scale2range,arParamRange);
    else
        p=randi(6)-1;
        x=generatePattern(p,T,range,scale2range);
    end
    % the whole series is already filled. just change it at appropriate
    % times
    if(~isempty(chgLocs))
        last=numel(chgLocs);
        for j=1:last
            l=chgLocs(j);
            chgs(l)=1;
            if(rand(1,1)<pAR)
                param=(arModelRange(2)-arModelRange(1))*randi(1,1)+arModelRange(1);
                tmp=generateRandAR(T-l+1,param,range,scale2range,arParamRange);      
                x(l:T)=tmp(1:T-l+1);                
            else
                p=randi(6)-1;
                tmp=generatePattern(p,T,range,scale2range);
                x(l:T)=tmp(1:T-l+1);                
            end
        end            
    end    

    if(~adjustEPA)
        scale2range=ss2range;
        range=srange;   
        x=adjust2range(x,range,scale2range);
    end
    
    if(nargout>3)
        xNoiseless=x;
    end
    
    if(smooth)
        x=smooth(x,ceil(0.2*rand(1,1)*minChangeDistance)*2+1);
    end
    
    if(noiseSigma>eps)
        x=x+noiseSigma.*randn(size(x));
    end
    if(noiseLevel>eps)
        x=x+((2*noiseLevel).*randn(size(x))-noiseLevel);
    end
    if pOutlier>eps
        if isempty(outlierRange)
            if scale2range
                if(sum(isinf(range))>0)
                    outlierRange=[-100,100];
                else
                    outlierRange=range;
                end
            else
                outlierRange=[-100,100];
            end
        end
        pout=rand(size(x));
        for i=1:numel(x)
            if(pout(i)<=pOutlier)
                x(i)=(outlierRange(2)-outlierRange(1))*rand(1,1)+outlierRange(1);
            end
        end
    end
end

function x=generatePattern(type,length,range,scale2range,param)
    x=zeros(length,1);
    length=length-1;    
    if(nargin<5)
        param=20*rand(1,1);
    end
    amplitude=10*(rand(1,1)-.5);
    rangeEnd=10*(rand(1,1)-.5);
    switch(type)
        case 0
            x=zeros(1,length+1);        
        case 1
            n=param;
            x=amplitude*sin(0:2*pi*n/double(length):2*pi*double(n));        
        case 2
            s=amplitude; e=rangeEnd;
            x=s:(e-s)/(length):e;
        case 3
            x=amplitude*ones(1,length+1);
        case 4
            n=length/param;
            p=double(param);
            x=amplitude*atan(0:2*pi/p:2*pi*double(n));        
        case 5
            n=length/param;
            x=amplitude*cos(0:2*pi/double(param):2*pi*double(n));        
        otherwise
            error('error : unknown pattern %d',type);
        %case 6
         %   x=generateRandAR(length,param);
    end
    x=adjust2range(x,range,scale2range);
end

function x=generateRandAR(length,n,range,scale2range,paramRange)
    d=poly((paramRange(2)-paramRange(1))*rand(1,n)+paramRange(1)); 
    x=generateAR(d,length,range,scale2range);
end

function x=generateAR(d,length,range,scale2range)
    d(2:end)=-d(2:end); 
    x=zeros(length,1); 
    x(1:size(d,2))=randn(size(d,2),1);
    for i=size(d,2)+1:size(x,1); 
        x(i)=d(2:end)*x(i-1:-1:i-size(d,2)+1); 
    end;
    x=adjust2range(x,range,scale2range);
end

function x=adjust2range(x,range,scale2range)
    if(isinf(range(1)))
        range(1)=min(x);
    end
    if(isinf(range(2)))
        range(2)=max(x);
    end    
    if(abs(max(x)-min(x))<eps)
        if(x(1)<range(1))
            x=range(1).*ones(size(x));
        end
        if(x(1)>range(2))
            x=range(2).*ones(size(x));
        end
    else
        if (~scale2range) && (max(x)<=range(2) && min(x)>=range(2))
           return;         
        end
        if(~scale2range) 
            if max(x)<=range(2); range(2)=max(x); end;
            if min(x)>=range(1); range(1)=min(x); end;
        end
        a=(range(1)-range(2))/(min(x)-max(x));
        b=range(1)-a*min(x);
        x=a.*x+b*ones(size(x));
    end
end
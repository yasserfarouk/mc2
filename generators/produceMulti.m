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

function [x,nDim,chgLocs,chgs]=produceMulti(M,T,varargin)
% prodcues example signals with patterns that are activated based on a
% causal model between n signals. A change in one signal will causally
% affect its cause signals.
%
% M             the model used as a cell of n*4 elements. The first is K*5
%               elements representing the parents of M in the causality
%               graph and the second is a single number
%               if M{i,1}(k,1)=j this means that change in i
%               is activated by a change in j (j-->i)
%               M{i,1}(k,2)=mean of delay in change in i due to change in j
%               M{i,1}(k,3)=std.dev of delay in change in i due to change in j
%               M{i,1}(k,4)=the probability of the causal connection
%               M{i,1}(k,5)=coupling number of parameters for the AR model.
%               this means the number of parameters of j in the new ar
%               model for i after the change. if zero then no coupling and
%               if negative then the new pattern in i need not be from an
%               AR model at all.
%               M{i,2}=the probability that i may change its dynamics
%               without any change in anything. 
%               M{i,3}= the range of AR parameters for each process
%               alone when the pattern is AR. 
%               M{i,4}= the number of dimensions for this time series. 
%
% T             The length of the time series
%
% pchangesNoParents  The probability that a change can happen in a process
%                    if it has no parent
%
%
%
% Outputs:
% ========
% x             The matrix  representing the time series (sum(nDim)*T)
% nDim          n*1 matrix representing the number of dimensions for each
%               signal.
%
%
    pAR=0.5;            % probability that the geneation is from an AR Model
    minPatternLength=ceil(T/20);
    %pchangeNoReason=0.0;
    pchangeFreeProcesses=0.01;
    range=[-inf, inf];
    scale2range=0;
    adjustEPA=1;
    nArgs=size(varargin,2);
    smooth=0;
    if(nArgs>0)
        if(mod(nArgs,2)~=0)
            error('The optional arguments must be in the form name,value so they must be even!!!');
        end
        for i=1:2:nArgs
            switch(lower(varargin{i}))       
                case {'pchangenoparent','pcnp','pchangefree'}        % probability of a change with absolutely no reason for processes that has no parent. default 0.01
                    pchangeFreeProcesses=(varargin{i+1});
                %case {'pchangenoparent','pcnp'}        % probability of a change with absolutely no reason. default 0
                 %   pchangeNoReason=(varargin{i+1});               
                case {'minpatternlength'}              % the minimum distance between any two allowable changes. default T/100. Notice that sometimes changes may happen faster than that if a time series had multiple parents (causes)
                    minPatternLength=(varargin{i+1});                
                case {'pautoregressive','par'}         % the probability of using an autoregressive model
                    pAR=(varargin{i+1});     
                case {'range'}                         % the range to allow the output to be within. Default -Inf:Inf
                    range=(varargin{i+1});     
                case {'scale2range'}                   % if nonzero then the output will be scaled to take the whole range form range(1) to range (2). default: 0
                    scale2range=(varargin{i+1});     
                case {'adjusteverypatternalone'}       % if nonzero each pattern will be adjusted to the range alone otherwise, the whole signal will be adjusted in the end. default:1
                    adjustEPA=(varargin{i+1});     
                case {'smooth'}                        % if nonzero the output will be smoothed with minpatternlength kernel. the default is no smoothing
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
    
    n=size(M,1);    
    nDim=ones(n,1);
    dRange=zeros(n,2);
    for i=1:n
        nDim(i)=M{i,4};
        if(i>1)
            dRange(i,1)=dRange(i-1,2)+1;
        else
            dRange(i,1)=1;
        end
        dRange(i,2)=dRange(i,1)+nDim(i)-1;
    end
    nt=sum(nDim);       % total dimensions
    x=zeros(nt,T);
    chgs=zeros(nt,T);
    chgLocs=cell(nt,1);
    done=boolean(zeros(n,1));
    nDone=0;
    % fill all time series that depend on nothing
    %j=1;
    
    % implement coupling
    for i=1:n
        for k=dRange(i,1):dRange(i,2)            
            % if self changes are allowed then put locations of self
            % change
            if(M{i,2}>eps)
                pchg=rand(1,T);
                l=1;
                while l<=T-minPatternLength
                    if(pchg(l)<M{i,2})                 
                        chgLocs{i,1}=[chgLocs{i,1};l];                        
                        l=l+minPatternLength-1;
                    end
                    l=l+1;
                end            
            end
        end
        % if that is all indicate that this dimension is done
        if(isempty(M{i,1}))            
            done(i)=1;
            nDone=nDone+1;        
            if(pchangeFreeProcesses>0)
                pcutoff=pchangeFreeProcesses-M{i,2};
                if(pcutoff>eps)
                    pchg=rand(1,T);
                    l=1;
                    while l<=T-minPatternLength
                        if(pchg(l)<pcutoff)                            
                            chgLocs{i,1}=[chgLocs{i,1};l];
                            l=l+minPatternLength-1;
                        end
                        l=l+1;
                    end            
                end
            end
        end
        %j=i+nDim(i);
    end
    % fill dependent series
    %trial=0;
    while nDone<n
        %trial=trial+1;
        %disp(sprintf('nDone:%02d  Trial:%03d\n',nDone,trial));
        for i=1:n
            if(~done(i))
                nI=size(M{i,1},1);
                mineDone=1;
                for j=1:nI
                    parent=M{i,1}(j,1);
                    if(done(parent)==0)
                        mineDone=0;
                    else
                        nch=size(chgLocs{parent,1});
                        for k=1:nch
                            if rand(1,1)<=M{i,1}(j,4)
                                tch=int32(round(chgLocs{parent,1}(k)+M{i,1}(j,3)*randn(1,1)+M{i,1}(j,2)));
                                if(tch>T-minPatternLength)
                                    tch=T-minPatternLength;
                                end
                                if(tch<1)
                                    tch=1;
                                end
                                chgLocs{i,1}=[chgLocs{i,1}; double(tch)];
                            end
                        end
                    end
                end
                if(mineDone)
                    done(i)=1;
                    nDone=nDone+1;
                end
            end
        end
    end    
    
    % do the changes in the correct time steps in all dimensions
    for i=1:nt        
        % fill with fixed pattern
        if(rand(1,1)<pAR)
            param=round((M{i,3}(2)-M{i,3}(1))*rand(1,1)+M{i,3}(1));
            tmp=generateRandAR(T,param,range,scale2range);                        
        else
            p=(round((5.49+0.49)*rand(1,1)-0.49));
            tmp=generatePattern(p,T,range,scale2range);
        end                    
        x(i,1:T)=tmp(1:T);            
        % the whole series is already filled. just change it at appropriate
        % times
        if(~isempty(chgLocs{i,1}))
            last=numel(chgLocs{i,1});
            for j=1:last
                l=chgLocs{i,1}(j);
                chgs(i,l)=1;
                if(rand(1,1)<pAR)
                    param=round((M{i,3}(2)-M{i,3}(1))*rand(1,1)+M{i,3}(1));
                    tmp=generateRandAR(T-l+1,param,range,scale2range);      
                    x(i,l:T)=tmp(1:T-l+1);                
                else
                    p=(round((5.49+0.49)*rand(1,1)-0.49));
                    tmp=generatePattern(p,T,range,scale2range);
                    x(i,l:T)=tmp(1:T-l+1);                
                end                    
                
            end            
        end
    end    
    
    if(~adjustEPA)
        scale2range=ss2range;
        range=srange;   
        for i=1:size(x,1)
            x(i,:)=adjust2range(x(i,:),range,scale2range);
        end
    end
    if(smooth)
        for i=1:size(x,1)
%            x(i,:)=smooth(x(i,:),ceil(0.2*rand(1,1)*minPatternLength)*2+1)
%            ; %todo make smooth work
        end
    end
end

function x=generatePattern(type,length,range,scale2range,param)
    x=zeros(length,1);
    length=length-1;    
    if(nargin<5)
        param=20*rand(1,1);
    end    
    switch(type)
        case 0
            x=zeros(1,length+1);        
        case 1
            n=param;
            x=sin(0:2*pi*n/double(length):2*pi*double(n));        
        case 2
            s=rand(1,1); e=rand(1,1);
            x=s:(e-s)/(length):e;
        case 3
            x=rand(1,1)*ones(1,length+1);
        case 4
            n=length/param;
            p=double(param);
            x=atan(0:2*pi/p:2*pi*double(n));
        case 5
            x=1:-2/double(length):-1; 
        case 6
            n=length/param;
            x=cos(0:2*pi/double(param):2*pi*double(n));        
        otherwise
            error('error : unknown pattern %d',type);
        %case 6
         %   x=generateRandAR(length,param);
    end
    x=adjust2range(x,range,scale2range);
end

function x=generateRandAR(length,n,range,scale2range)
    d=poly(2*rand(1,n)-1); 
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
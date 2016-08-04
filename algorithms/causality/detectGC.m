function [M,F,cvUsed]=detectGC(x,alpha,max_lag,varargin)
% applies granger_cause (granger causality test with BIC for optimal lag
% selection) to all the dimensions of X and generates the causality model
% from which this data was generated
%
%
% Inputs:
% ======
% x             n*T array of time series to find the relations between its
%               n timeseries
% alpha         The significance level to be used to decide causality
% max_lag       the maximum lag used in testing causality (BIC will be used
%               to select the optimal causality less than or equal to this)
% varargin      see the switch statement
% 
% Outputs:
% =======
% M             The model (look at produceMulti to understand its
%               structure)
% F             n*n array of The F values returned by the granger causality
%               test (i,j) is the result of testing that j g-causes i
% cvUsed        n*n array of the critical values with which F was compared.
%               This is by default is coming from F distribution but if
%               'th' was specified in the varargin, it will be th. IF you
%               supply a vector threshold then it will be filled with it

findSelfLoops=boolean(0);
n=size(x,1);
th=-1;
if(size(x,2)<n)
    warning ('the number of time steps is less than the number of dimension... confirm that you did not forget a transpose');
end

nArgs=size(varargin,2);
if(nArgs>0)
    if(mod(nArgs,2)~=0)
        error('The optional arguments must be in the form name,value so they must be even!!!');
    end
    for i=1:2:nArgs
        switch(lower(varargin{i}))       
            case {'findselfloops','fsl'} 
                findSelfLoops=(varargin{i+1});            
            case {'th','threshold'}             % if not negative then rather than comparing F with the critical value from Fstatistic we compare it with th.
                th=(varargin{i+1});             % if th is a vector of m elements then Mg will be an array of m elements containing all the models with different thresholds
            otherwise
                error('Unknown argument: %s',varargin{i});
        end
    end
end

F=zeros(n,n);
if(th(1)<0)    
    m=numel(alpha);
    cvUsed=zeros(n,n,m);
    M=cell(n,m);      
    for i=1:n
        for j=1:n
            if (findSelfLoops==0) && (i==j)
                continue;
            end
            calculated=0;
            for k=1:m        
                if(alpha(k)<eps)
                    continue;
                end
                if(~calculated)
                    calculated=1;
                    [F(i,j),cvUsed(i,j,k),x_lag,y_lag,T]=granger_cause(x(i,:),x(j,:),alpha(k),max_lag);                    
                end
                if(alpha(k)>=1)
                    M{i,k}= [M{i,k}; j,y_lag,0,1,y_lag];
                    continue;
                end
                cvUsed(i,j,k) = finv(1-alpha(k),y_lag,(T-(x_lag+y_lag+1)));                
                if(F(i,j)>cvUsed(i,j,k))                    
                    M{i,k}= [M{i,k}; j,y_lag,0,1,y_lag];
                end            
            end
        end
    end                
else
    m=numel(th);
    cvUsed=zeros(n,n,m);
    M=cell(n,m);    
    for i=1:n
        for j=1:n
            if (findSelfLoops==0) && (i==j)
                continue;
            end
            [F(i,j),cv,x_lag,y_lag]=granger_cause(x(i,:),x(j,:),0.05,max_lag);                            
            for k=1:m
                cvUsed(i,j,k)=th(k);
                if(F(i,j)>th(k))                    
                    M{i,k}= [M{i,k}; j,y_lag,0,1,y_lag];
                end        
            end
        end
    end            
end
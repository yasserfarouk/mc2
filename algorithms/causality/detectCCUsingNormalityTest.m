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

function [M,P]=detectCCUsingNormalityTest(locs,th,varargin)
% finds causality structure using normality of delays test
%
% Inputs:
% ======
% locs          a n*1 cell array each of them containing the change locations
%               associated with one time series
% th            the threshold used to get the granger causality. can be set
%               to a vector of length m and in this case the output M will be
%               a cell array of n*m members where M(:,i) corresponds to
%               th(i)
% varargin      see the swich statement
%
%
%
%
%
%
%
%
%
%
findSelfLoops=boolean(0);
nDelayTrials=1;
rmBottom=0.1;
rmTop=0.1;
maxLag=bitmax;
nArgs=size(varargin,2);
if(nArgs>0)
    if(mod(nArgs,2)~=0)
        error('The optional arguments must be in the form name,value so they must be even!!!');
    end
    for i=1:2:nArgs
        switch(lower(varargin{i}))       
            case {'findselfloops','fsl'} 
                findSelfLoops=(varargin{i+1});  
            case {'maxlag','maxdelay'}              % The maximum delay to look for relations in the time series dynamics (default 10)
                maxLag=(varargin{i+1});                        
            case {'ndelaytrials','nd'}          % the number of time delays to calculate for every change point. default is bitmax (disabled)
                nDelayTrials=(varargin{i+1});    
            case {'rmbottom','rb'}              % the fraction of time delays to ignore in the bottom of the list. default is 0
                rmBottom=(varargin{i+1});       
            case {'rmtop','rt'}                 % the fraction of time delays to ignore in the top of the list. default is 0
                rmTop=(varargin{i+1});       
            otherwise
                error('Unknown argument: %s',varargin{i});
        end
    end
end
n=length(locs);
m=numel(th);
M=cell(n,m);
P=zeros(n,n,m);
tau=cell(n,1);
for i=1:n    
    if(findSelfLoops)
        nt=numel(locs{i});
        if(nt>0)
            for j=1:nt-1
                ndt=min(nDelayTrials,nt-j);
                for k=j+1:j+ndt
                  tau{i}= [tau{i};locs{i}(k)-locs{i}(j)];
                end
            end
        end
    end
end
for k=1:m   
    if(th(k)<eps)
        continue;
    end
    for i=1:n
        for j=1:n
            if (i==j)
                if findSelfLoops
                    ds=tau{i};
                    nt=numel(ds);
                    if(nt>0)
                        nr=round(rmBottom*nt);
                        if(nr>0)
                            ds(1:nr)=[];            
                        end
                        nr=round(rmTop*nt);
                        if(nr>0)
                            ds(nt-nr+1:end)=[];            
                        end        
                    end
                    ds(ds>maxLag)=[];
                    nt=numel(ds);
                else
                    nt=0;            
                end
            else
                ds=getDelays(locs{j},locs{i},nDelayTrials,rmBottom,rmTop);
                ds(ds>maxLag)=[];
                nt=numel(ds);
            end            
            if(th(k)>=1 && ~(i==j && ~findSelfLoops))
              if isempty(ds); 
                  mds=NaN; vars=Inf; 
              else                  
                  mds=mean(ds); vars=std(ds);
              end;
              M{i,k}= [M{i,k}; j,mds,vars,0,0];              
              continue;
            end
            if(nt>3)
                ds=double(sort(ds));
                [DP,P(i,j,k)] = DagosPtest(ds,th(k));
                if(DP)                
                    M{i,k}= [M{i,k}; j,mean(ds),std(ds),1-P(i,j,k),0];
                end     
            end
        end
    end
end
end

function ds=getDelays(x,y,ndt,rmBottom,rmTop)
% calculates the delays from changes in x to changes in y
    ds=[];
    n=numel(x); m=numel(y);
    for i=1:n
        fnd=0;
        for j=1:m
            if(y(j)>=x(i))
                ds=[ds;y(j)-x(i)];
                fnd=fnd+1;
                if(fnd>=ndt)
                    break;
                end
            end
        end
    end
    n=numel(ds);
    if(n>0)
        nr=round(rmBottom*n);
        if(nr>0)
            ds(1:nr)=[];            
        end
        nr=round(rmTop*n);
        if(nr>0)
            ds(n-nr+1:end)=[];            
        end        
        ds=double(ds);
    end
end
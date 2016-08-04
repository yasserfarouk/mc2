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

function [locs,locsAboveThreshold]=findLocsTh(r,varargin)
% first parameter is the time series in which we want to localize local
% maxima
% optional parameters in order
% th        The threshold underwhich we do not consider any thing (default
%           is .05*max(r)). If multiple thresholds are given then the result is
%           a cell array of changes with each threshold alone
% firstW    if firstW<-1 then the location of the local maxima is just returned
%           if 0>firstW>=-1 then one of the points to attain the local maxima is
%           to be used as the location of the maxima using:
%           firstW*firstLocalMaxima+(1-firstW)*lastLocalMaxima
%           where firstLocalMaxima(lastLocalMaxima) are the location within
%           a local vecinity of a local maxima that are within th value
%           from it
%           if positive then the location of the maxima is calculated as:
%           firstW*first+(1-firstW)*last 
%           where first and last are the points at which the signal goes
%           over and returns under the threshold th
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


    th=.05*max(r);
    firstW=-2;     % default is location of maximum change
    nArgs=size(varargin,2);
    if(nArgs>0)
        if ~isempty(varargin{1})
            th=varargin{1};
        end
    end
    if(nArgs>1)
        if ~isempty(varargin{1})
            firstW=varargin{2};        
        end
    end
    JustMax=0;
    if(firstW<-1-eps)
        firstW=0;
        lastW=0;
        useMaxPt=1;
        JustMax=1;
    elseif(firstW>=-eps)    
        lastW=1-firstW;
        useMaxPt=0;
    else
        firstW=-firstW;
        lastW=1-firstW;
        useMaxPt=1;
    end
    m=length(th);
    if m==1
        [locs,locsAboveThreshold]=doFindLocs(r,th,firstW,lastW,useMaxPt,JustMax);
    else
        locs=cell(1,m);
        locsAboveThreshold=cell(1,m);
        for t=1:m
            [locs{t},locsAboveThreshold]=doFindLocs(r,th(t),firstW,lastW,useMaxPt,JustMax);
        end
    end    
end

function [locs,locsAboveThreshold]= doFindLocs(r,th,firstW,lastW,useMaxPt,JustMax)
    K=numel(r);    
    locs=[];
    in1=boolean(0);
    first1=0;    
    maxPt=0;
    firstMax=0;
    lastMax=0;
    inMax=boolean(0);
    locsAboveThreshold=[];
    for i=1:K
        if in1
            if r(i)>r(maxPt)
                maxPt=i;
            end
            if r(i)<th
                last1=i-1;
                in1=0;
                locsAboveThreshold(end,2)=i-1;
                if useMaxPt
                    if JustMax
                        locs=[locs;maxPt];
                    else
                        maxMin=r(maxPt)-th;
                        for j=first1:last1
                            if r(j)>=maxMin
                                if inMax
                                    lastMax=j;
                                else
                                    inMax=1;
                                    firstMax=j;
                                    lastMax=firstMax;
                                end
                            end
                        end
                        inMax=0;
                        locs=[locs;round(firstW*firstMax+lastW*lastMax)];
                    end
                else                    
                    locs=[locs;round(firstW*first1+lastW*last1)];
                end                
            end            
        else
            if r(i)>=th                            
                first1=i;
                last1=first1;
                in1=1;
                maxPt=i;
                locsAboveThreshold=[locsAboveThreshold;i,-1];
            end
        end        
    end    
    if in1
        in1=0;
        last1=i;
        if useMaxPt
            if JustMax
                locs=[locs;maxPt];
            else
                maxMin=r(maxPt)-th;
                for j=first1:last1
                    if r(j)>=maxMin
                        if inMax
                            lastMax=j;
                        else
                            inMax=1;
                            firstMax=j;
                            lastMax=firstMax;
                        end
                    end
                end
                inMax=0;
                locs=[locs;round(firstW*firstMax+lastW*lastMax)];
            end
        else
            locs=[locs;round(firstW*first1+lastW*last1)];
        end
    end
    if ~isempty(locsAboveThreshold)
        if (locsAboveThreshold(end,2)<0)
            locsAboveThreshold(end,2)=K;
        end
    end
end

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

function [begends,toreuse]=extendStem(nPnts,x,T,maxLen,begends,...    
    stopFun,toreuse,stopFunParams,...
    maxSmallDist,minLargeDist,dFun,dFunParams,...
    minExtension,tryExtendingStemsSlowly,useBisection...
    ,singleTrial,ignoreOverlapsDuringExtension)
% Extends a motif stem (internal function)
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%

    stop=0;    
    if nargin<16
        singleTrial=0;
    end
    nAdded=nPnts(2)-nPnts(1);
    if nAdded<minExtension
        if singleTrial
            return;
        else
            stop=1;            
        end
    end                    
    %newbegends=[max([ones(size(begends,1),1),begends(:,1)+nPnts(1)],[],2),...
    %    min([T.*ones(size(begends,1),1),begends(:,2)+nPnts(2)],[],2)];
    newbegends=[begends(:,1)+nPnts(1),begends(:,2)+nPnts(2)];
    if (min(newbegends(:,1))<1) || (max(newbegends(:,2))>T)
        if singleTrial
            return;
        else
            stop=1;            
        end
    end
    if max(max(abs(newbegends-begends)))<1;
        if singleTrial
            return;
        else
            stop=1;            
        end
    end
    if max(newbegends(:,2)-newbegends(:,1)+1)>maxLen
        stop=1;        
    end
    if ~stop && ~isinf(minLargeDist)
        stop=isMeanDistBetweenOccurencesSmall(x,newbegends,maxSmallDist,minLargeDist,dFun,dFunParams);
    end
    if(~stop)
        [stop,toreuse]=stopFun(x,begends,newbegends,toreuse,stopFunParams);
        if ~stop; begends=newbegends; end;
    end
    if tryExtendingStemsSlowly
        while(~stop)               
            begends=newbegends;
            nStemsInGroup=size(begends,1);            
            % check that no occurrences overlap            
            if ~ignoreOverlapsDuringExtension
                overlap=[];
                for j=1:nStemsInGroup
                    for k=1:nStemsInGroup
                        if(k==j) 
                            continue;
                        end
                        if testOverlapping(begends(j,:),begends(k,:))
                            overlap=[overlap;k,j];
                        end
                    end
                end
                % if there is an overlap, we just stop and exit
                if ~isempty(overlap)
                    % @todo do something better in case of overlap!!!!!!!
                    %warning('DGR:OVERPLAP','Overlap detected ... DGR');
                    break;
                end
            end

            %newbegends=[max([ones(size(begends,1),1),begends(:,1)+nPnts(1)],[],2),...
            %    min([T.*ones(size(begends,1),1),begends(:,2)+nPnts(2)],[],2)];
            newbegends=[begends(:,1)+nPnts(1),begends(:,2)+nPnts(2)];
            if (min(newbegends(:,1))<1) || (max(newbegends(:,2))>T)
                if singleTrial
                    return;
                else
                    stop=1;            
                end
            end
            if max(max(abs(newbegends-begends)))<1;
                if singleTrial
                    return;
                else
                    stop=1;            
                end
            end
            if max(newbegends(:,2)-newbegends(:,1)+1)>maxLen
                if singleTrial
                    return;
                else
                    stop=1;            
                end
            end
            if ~stop && ~isinf(minLargeDist)
                stop=isMeanDistBetweenOccurencesSmall(x,newbegends,maxSmallDist...
                    ,minLargeDist,dFun,dFunParams);
            end
            if(~stop)
                [stop,toreuse]=stopFun(x,[],newbegends,toreuse,stopFunParams);
            end
        end
    end
    % this means last trial failed. we will try to use smaller nPnts at
    % each side until it is not possible anymore
    if stop
        if useBisection
            nPnts2Try=[nPnts(1),0;0,nPnts(2)];
            for d=1:2            
                nPntsMax=double(nPnts2Try(d,:));
                nPntsMin=[0,0];
                while(1)
                    stop=0;
                    if max(abs(nPntsMax-nPntsMin))<=minExtension
                        if any(nPntsMin)
                            %begends=[max([ones(size(begends,1),1),begends(:,1)+nPntsMin(1)],[],2),...
                            %    min([T.*ones(size(begends,1),1),begends(:,2)+nPntsMin(2)],[],2)];
                            newbegends=[begends(:,1)+nPntsMin(1),begends(:,2)+nPntsMin(2)];
                            if ~((min(newbegends(:,1))<1) || (max(newbegends(:,2))>T))
                               begends=newbegends;
                            end
                        end
                        break;
                    end                
                    nPnts=round((nPntsMax+nPntsMin)./2);
                    nAdded=nPnts(2)-nPnts(1);
                    if nAdded<minExtension
                        if any(nPntsMin)
                            %begends=[max([ones(size(begends,1),1),begends(:,1)+nPntsMin(1)],[],2),...
                            %    min([T.*ones(size(begends,1),1),begends(:,2)+nPntsMin(2)],[],2)];
                            newbegends=[begends(:,1)+nPntsMin(1),begends(:,2)+nPntsMin(2)];
                            if ~((min(newbegends(:,1))<1) || (max(newbegends(:,2))>T))
                               begends=newbegends;
                            end
                        end
                        break;
                    end                
                    %newbegends=[max([ones(size(begends,1),1),begends(:,1)+nPnts(1)],[],2),...
                    %    min([T.*ones(size(begends,1),1),begends(:,2)+nPnts(2)],[],2)];
                    newbegends=[begends(:,1)+nPnts(1),begends(:,2)+nPnts(2)];
                    if (min(newbegends(:,1))<1 || max(newbegends(:,2))>T)
                       stop=1;
                    end
                    if max(max(abs(newbegends-begends)))<1;
                        stop=1;
                    end
                    if max(newbegends(:,2)-newbegends(:,1)+1)>maxLen
                        stop=1;        
                    end
                    if ~stop && ~isinf(minLargeDist)
                        stop=isMeanDistBetweenOccurencesSmall(x,newbegends,maxSmallDist,minLargeDist,dFun,dFunParams);
                    end
                    if(~stop)
                        [stop,toreuse]=stopFun(x,begends,newbegends,toreuse,stopFunParams);
                    end
                    if stop
                        nPntsMax=nPnts;
                    else
                        nPntsMin=nPnts;                    
                    end                
                end
            end
        else
            if singleTrial
                return;
            end
            if nPnts(1) ~=0 && nPnts(2) ~=0
                [begends,toreuse]=tryExtendingBefore(nPnts,x,T,maxLen,begends,...    
                    stopFun,toreuse,stopFunParams,...
                    maxSmallDist,minLargeDist,dFun,dFunParams,...
                    minExtension,tryExtendingStemsSlowly,useBisection,ignoreOverlapsDuringExtension);                
                [begends,toreuse]=tryExtendingAfter(nPnts,x,T,maxLen,begends,...    
                    stopFun,toreuse,stopFunParams,...
                    maxSmallDist,minLargeDist,dFun,dFunParams,...
                    minExtension,tryExtendingStemsSlowly,useBisection,ignoreOverlapsDuringExtension);
            elseif nPnts(1) ~=0
                [begends,toreuse]=tryExtendingBefore(nPnts,x,T,maxLen,begends,...    
                    stopFun,toreuse,stopFunParams,...
                    maxSmallDist,minLargeDist,dFun,dFunParams,...
                    minExtension,tryExtendingStemsSlowly,useBisection,ignoreOverlapsDuringExtension);
            else
                [begends,toreuse]=tryExtendingAfter(nPnts,x,T,maxLen,begends,...    
                    stopFun,toreuse,stopFunParams,...
                    maxSmallDist,minLargeDist,dFun,dFunParams,...
                    minExtension,tryExtendingStemsSlowly,useBisection,ignoreOverlapsDuringExtension);
            end
        end
    end
end
function [begends,toreuse]=tryExtendingBefore(nPnts,x,T,maxLen,begends,...    
    stopFun,toreuse,stopFunParams,...
    maxSmallDist,minLargeDist,dFun,dFunParams,...
    minExtension,tryExtendingStemsSlowly,useBisection,ignoreOverlapsDuringExtension)
        [begends,toreuse]=extendStem([nPnts(1),0],x,T,maxLen,begends,...    
            stopFun,toreuse,stopFunParams,...
            maxSmallDist,minLargeDist,dFun,dFunParams,...
            minExtension,tryExtendingStemsSlowly,useBisection,1,ignoreOverlapsDuringExtension);
        [begends,toreuse]=extendStem([sign(nPnts(1))*minExtension,0],x,T,maxLen,begends,...    
            stopFun,toreuse,stopFunParams,...
            maxSmallDist,minLargeDist,dFun,dFunParams,...
            minExtension,tryExtendingStemsSlowly,useBisection,1,ignoreOverlapsDuringExtension);
end
function [begends,toreuse]=tryExtendingAfter(nPnts,x,T,maxLen,begends,...    
    stopFun,toreuse,stopFunParams,...
    maxSmallDist,minLargeDist,dFun,dFunParams,...
    minExtension,tryExtendingStemsSlowly,useBisection,ignoreOverlapsDuringExtension)
        [begends,toreuse]=extendStem([nPnts(1),0],x,T,maxLen,begends,...    
            stopFun,toreuse,stopFunParams,...
            maxSmallDist,minLargeDist,dFun,dFunParams,...
            minExtension,tryExtendingStemsSlowly,useBisection,1,ignoreOverlapsDuringExtension);
        [begends,toreuse]=extendStem([sign(nPnts(1))*minExtension,0],x,T,maxLen,begends,...    
            stopFun,toreuse,stopFunParams,...
            maxSmallDist,minLargeDist,dFun,dFunParams,...
            minExtension,tryExtendingStemsSlowly,useBisection,1,ignoreOverlapsDuringExtension);
end
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

function [locs,means,stats,maxSmallDistance,execTime]=mknRange(x,lrange,lStep,nMotifs,Fast)
% candLocs is not used
   % lStep=min((lrange(2)-lrange(1)),max(4,ceil((lrange(2)-lrange(1))/10)));       
   locs=cell(0,0);
   means=cell(0,0);
   stats=cell(0,0);
   execTime=0;
   for L=lrange(1):lStep:lrange(2)
        [locst,exect,~,meanst,statst,maxSmallDistance]=mkn(x,L,0.0,10,nMotifs,0.0,0,Fast);%,lStep,innerOverlap,R ...
                        %,nMotifs,outerOverlap, Fast,recalcDists);
        execTime=execTime+exect;
        if(size(locst,1)>0)
            %[locst,meanst,statst]=mkCombine(x,locst,L);
            locs=[locs;locst];
            means=[means,meanst];
            stats=[stats,statst];
        end
   end
   
end
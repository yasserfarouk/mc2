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

function [locs,means,stats,maxSmallDistance]=mknnC(x,lrange,candLocs)
%candLocs is not used
   % lStep=min((lrange(2)-lrange(1)),max(4,ceil((lrange(2)-lrange(1))/10)));    
   
   [locs,dists]=mknn(x,lrange,1);%,lStep,innerOverlap,R ...
   n=length(locs);
   L=length(x)-lrange(2);
   toremove=[];
   for i=1:n
       if locs{i}(1,2)>L || locs{i}(2,2)>L
           toremove=[toremove,i];
       end
   end
   if ~isempty(toremove)
       locs(toremove)=[];
       dists(toremove)=[];
       means(toremove)=[];
       stats(toremove)=[];
   end
                        %,nMotifs,outerOverlap, Fast,recalcDists);
   [locs,means,stats,~,maxSmallDistance]=mkCombinePerLength(x,locs,dists);
   
end
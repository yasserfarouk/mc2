function [locs,means,stats,maxSmallDistance]=mnC(x,lrange,candLocs,nMotifs)
%candLocs is not used
   % lStep=min((lrange(2)-lrange(1)),max(4,ceil((lrange(2)-lrange(1))/10)));    
   if nargin<4
       nMotifs=10;
   end
   [locs,dists]=mn(x,lrange,5,3,0,nMotifs);%,lStep,innerOverlap,R ...
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
   end
                        %,nMotifs,outerOverlap, Fast,recalcDists);
   [locs,means,stats,~,maxSmallDistance]=mkCombinePerLength(x,locs,dists);
   
end
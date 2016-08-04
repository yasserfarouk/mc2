function [locs,means,stats,maxSmallDistance]=mkpC(x,lrange,candLocs)
% candLocs is not used
   % lStep=min((lrange(2)-lrange(1)),max(4,ceil((lrange(2)-lrange(1))/10)));    
   lStep=max(10,ceil((lrange(2)-lrange(1))/10));
   locs=cell(0,0);
   means=cell(0,0);
   stats=cell(0,0);
   dists=[];
   for L=lrange(1):lStep:lrange(2)
        [locst,ds,~,meanst,statst]=mkp(x,L,0.0,8,15,0.0,0,0);%,lStep,innerOverlap,R ...
                        %,nMotifs,outerOverlap, Fast,recalcDists);
        if(size(locst,1)>0)
            %[locst,meanst,statst]=mkCombine(x,locst,L);
            locs=[locs;locst];
            means=[means,meanst];
            stats=[stats,statst];
            dists=[dists,ds];
        end
   end
   [locs,means,stats,~,maxSmallDistance]=mkCombinePerLength(x,locs,dists);
end
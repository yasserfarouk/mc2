function [locs,means,stats,maxSmallDistance,execTime]=mkpRange(x,lrange,lStep,nMotifs,Fast)
% candLocs is not used
   % lStep=min((lrange(2)-lrange(1)),max(4,ceil((lrange(2)-lrange(1))/10)));       
   locs=cell(0,0);
   means=cell(0,0);
   stats=cell(0,0);
   execTime=0;
   %strt=tic;
   for L=lrange(1):lStep:lrange(2)
        [locst,exect,~,meanst,statst,maxSmallDistance]=mkp(x,L,0.0,10,nMotifs,0.0,0,Fast);%,lStep,innerOverlap,R ...
                        %,nMotifs,outerOverlap, Fast,recalcDists);
        execTime=execTime+exect;
        if(size(locst,1)>0)
            %[locst,meanst,statst]=mkCombine(x,locst,L);
            locs=[locs;locst];
            means=[means,meanst];
            stats=[stats,statst];
        end
   end
   %execTime=1e6*toc(strt);
end
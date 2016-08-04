function [locs,means,stats,maxSmallDistance,execTime,errorLocs]=mkmmRange(x,lrange,lStep,...
    nMotifs,recalcDists)
% a wrapper to mkpp to make it look as if taking the same parameters as
% mkpRange ... used internally by mkTest
  
   [locs,execTime,~,means,stats,maxSmallDistance]=mkmm(x,lrange,lStep...
       ,0,10,nMotifs,0,recalcDists);%,lStep,innerOverlap,R ...
   n=length(locs);
   L=length(x)-lrange(2);
   toremove=[];
   for i=1:n
       if locs{i}(1,2)>L || locs{i}(2,2)>L
           toremove=[toremove,i];
       end
   end
   if ~isempty(toremove)
       errorLocs=locs(toremove);
       locs(toremove)=[];
       means(toremove)=[];
       stats(toremove)=[];
   end
                        %,nMotifs,outerOverlap, Fast,recalcDists);
   %[locs,means,stats]=mkCombine(x,locs,lrange(1));
   
end
function [locs,toremove]=mdRemoveShortMotifs(locs,Lmin)
n=numel(locs);
L=zeros(n,1);
for i=1:n
    if ~isempty(locs{i})
        L(i)=locs{i}(1,2)-locs{i}(1,1)+1;
    end
end
toremove=find(L<Lmin);
locs(toremove)=[];
end
function [locs,toremove]=mdRemoveNonMaximal(locs)
% removes nonmaximal motifs. Assumes that the motifs are stored in order of
% importance so it keeps the earliest possible motifs in this list
n=numel(locs);
toremove=[];
for i=1:n    
    for j=i+1:n
        if isincluded(locs{j},locs{i})
            toremove=[toremove j];
        end
    end
end

locs(toremove)=[];
end

function v=isincluded(loc1,loc2)
    n=size(loc1,1);
    m=size(loc2,1);
    v=zeros(n,1);
    for i=1:n
        for j=1:m
            if loc1(i,1)>=loc2(j,1) && loc1(i,2)<=loc2(j,2)
                v(i)=1;
            end
        end
    end
    v=all(v);
end
function [ locsF,meansF,statsF,dists,maxSmallDist ] = mkCombine(x, locs, L,dists)
%Combines the results of MK algorithm runs into sets of motifs. All
%returned motifs will be of length L.
%
% This function is intended to be used for combining motif pairs returned
% from exact motif discovery algorithms at some predefined length. It can
% be used though with pairs of different lengths but it will trim the ones
% longer than L when returning the results.


if isempty(locs)
    locsF=cell(0,1);
    meansF=cell(0,1);
    statsF=cell(0,1);
    maxSmallDist=0;
    return;
end
minTrivialOverlap=0.95;     % any overlap over this means that we are looking at the same thing
minIntersection=0.00;       % accept any intersection between motif sets
maxInternalOverlap=0.00;    % do not allow any internal overlap within occurrences in a set
if nargin<4
    dists=zeros(size(locs,1),1);
end
if iscell(locs)
    [maxSmallDist,locs,dists]=mkRemoveLargeDists(locs,dists);
    nPairs=size(locs,1);
    locsM=cell2mat(locs);
else
    maxSmallDist=max(dists);
    nPairs=size(locs,1)/2;
    locsM=locs;
end
nPoints=size(locsM,1);
% make the first column an occurrence beginning and the second its pair
for i=1:2:nPoints
    locsM(i,2)=i+1;
    locsM(i+1,2)=i;
end
overlappedWith=cell(nPoints,1);
%changeWith=zeros(nPoints,1);
pairedWith=cell(nPoints,1);
for i=1:nPoints
    pairedWith{i}=[pairedWith{i},locsM(i,2)];
end
begRange=[locsM(:,1),locsM(:,1)];
% replace each set of occurrences that are overlapping more than
% minTrivialOverlap with the occurrence in their center.
% TODO I can take the distances into account and shift to the lowest
% distances
somethingChanged=1;

while somethingChanged==1
    somethingChanged=0;
    replaced=zeros(nPoints,1);
    for i=1:nPoints-1
        %if replaced(i);      continue;    end
        for j=i+1:nPoints
            of=overlapfraction(...
                    [locsM(i,1),locsM(i,1)+L-1],...
                    [locsM(j,1),locsM(j,1)+L-1]...
                    ,1);
            if(of>=minTrivialOverlap)
                overlappedWith{i}=union(overlappedWith{i},[j]);
                pairedWith{i}=union(pairedWith{i},[locsM(j,2)]);
                begRange(i,1)=min(begRange(i,1),locsM(j,1));
                begRange(i,2)=max(begRange(i,2),locsM(j,1));            
                %if changeWith(j)==0
                %    changeWith(j)=i;
                %end
            end
        end
        nOverlapping=length(overlappedWith{i});
        if nOverlapping>0
            rng=begRange(i,2)-begRange(i,1);
            if rng>0                
                somethingChanged=1;
                rng=round(0.5*(begRange(i,1)+begRange(i,2)));            
                begRange(i,2)=rng;
                begRange(i,1)=rng;
                locsM(i,1)=rng;            
                for k=1:nOverlapping
                    j=overlappedWith{i}(k);
                    locsM(j,1)=rng;                
                    replaced(j)=1;
                end
            end
        end
    end
end
for i=1:length(pairedWith)
    n=length(pairedWith{i});
    for j=1:n
        pairedWith{i}(j)=locsM(pairedWith{i}(j),1);
    end
end
for i=1:2:nPoints
    pairedWith{i}=union(pairedWith{i},pairedWith{i+1});
end
pairedWith=pairedWith(1:2:nPoints);
%combine motifs that can be combined
tobermoved=zeros(nPairs,1);
for i=1:nPairs-1
    if tobermoved(i); continue; end;
    for j=i+1:nPairs
        if tobermoved(j); continue; end;
        inter=intersect(pairedWith{i},pairedWith{j});
        if ~isempty(inter)
            un=union(pairedWith{i},pairedWith{j});
            if length(inter)/length(un)>minIntersection;
                pairedWith{i}=un;
                tobermoved(j)=1;
            end
        end
    end
end

% remove motifs that are included in other motifs
nM=length(pairedWith);
for i=1:nM-1
    if tobermoved(i); continue; end;
    for j=i+1:nM
        if tobermoved(j); continue; end;
        if isempty(setxor(pairedWith{i},pairedWith{j}));
            pairedWith{i}=union(pairedWith{i},pairedWith{j});
            tobermoved(j)=1;
        end
    end
end
pairedWith(tobermoved==1)=[];
nM=length(pairedWith);
locsF=cell(nM,1);
meansF=[]; statsF=[];
for i=1:nM
    nO=length(pairedWith{i});
    for j=1:nO
        locsF{i}=[locsF{i};pairedWith{i}(j),...
            pairedWith{i}(j)+L-1];
    end
    [statsF,meansF]=updateStatsAndMean(locsF{i},x,@abs,statsF,meansF,i);
end

% order motifs to put the ones explaining more of the data first
[locsF,meansF,statsF]=mdOrder(locsF,meansF,statsF,-1);


% make occurrences in the same motif
torem=zeros(nM,1);
for i=1:nM-1
    nI=length(locsF{i});
    ovMat=zeros(nI,nI);
    begRange=[locsF{i}(:,1),locsF{i}(:,1)];
    overlappingWith=cell(nI,1);
    replaced=zeros(nI,1);
    for k=1:nI-1        
        for l=k+1:nI
            of=overlapfraction(...
                [locsF{i}(k,1),locsF{i}(k,2)],...
                [locsF{i}(l,1),locsF{i}(l,2)]...
                ,1);
            if of>maxInternalOverlap
                begRange(k,1)=min([begRange(k,1),locsF{i}(l,1)]);
                begRange(k,2)=max([begRange(k,1),locsF{i}(l,1)]);
                overlappingWith{k}=[overlappingWith{k},l];
            end
        end
        nOverlapping=length(overlappingWith{k});
        if nOverlapping>0
            rng=begRange(k,2)-begRange(k,1);
            if rng>0                
                somethingChanged=1;
                rng=round(0.5*(begRange(k,1)+begRange(k,2)));            
                begRange(k,2)=rng;
                begRange(k,1)=rng;
                locsF{i}(k,1)=rng;                        
                for l=1:nOverlapping
                    j=overlappingWith{k}(l);
                    locsF{i}(j,1)=rng;                                    
                    replaced(j)=1;
                end
            end
        end
    end
    tmp=unique(locsF{i}(:,1));
    locsF{i}=[tmp,tmp+L-1];
end

nM=length(locsF);
for i=1:nM    
    [statsF,meansF]=updateStatsAndMean(locsF{i},x,@abs,statsF,meansF,i);
end
[locsF,meansF,statsF]=mdOrder(locsF,meansF,statsF,-1);
end


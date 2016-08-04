function [locs,means,stats] = gemoda(x, Lmin, th,minSupport, distFun...
    ,isSymmetricDistFun,clusterFun, inLocs, minDistBetweenOccurenceBeg,Lmax)
%GEMODA applies GEMODA algorithm for motif discovery
%   


if isvector(x)
    x=x(:);
end
[T,n]=size(x);
set(0,'RecursionLimit',max([1000,T]));

if ~exist('minSupport','var') || isempty(minSupport)
    minSupport=2;
end
if ~exist('Lmax','var') || isempty(Lmax)
    Lmax=T;
else
    Lmax=Lmax+1;
end
if ~exist('clusterFun','var') || isempty(clusterFun)
    %clusterFun=@findCliques;
    clusterFun=@findConnectedComp;
end
if ~exist('distFun','var') || isempty(distFun)
    distFun=@distEuclidean;
    isSymmetricDistFun=true;
end
if ~exist('isSymmetricDistFun','var') || isempty(isSymmetricDistFun)
    isSymmetricDistFun=false;
end
if ~exist('minDistBetweenOccurenceBeg','var') || isempty(minDistBetweenOccurenceBeg)
    minDistBetweenOccurenceBeg=1;
end
if ~exist('inLocs','var') || isempty(inLocs)
    inLocs=1:T-Lmin+1;
else
    if numel(inLocs)==1
        inLocs=1:inLocs:T-Lmin+1;
    end
end
inLocs(inLocs<1)=[];
inLocs(inLocs>T-Lmin+1)=[];
nLocs=numel(inLocs);
% comparison
%rows=[1:T-Lmin+1]'; cols=[1:T-Lmin+1]';
rows=[]; cols=[];
for i=1:nLocs    
    for j=i+1:nLocs
        if (inLocs(j)-inLocs(i)<minDistBetweenOccurenceBeg)
            continue;
        end
        d=distFun(x(inLocs(i):inLocs(i)+Lmin-1,:),x(inLocs(j):inLocs(j)+Lmin-1,:));
        if d<=th
            rows=[rows;inLocs(i)];
            cols=[cols;inLocs(j)];
            if isSymmetricDistFun
                rows=[rows;inLocs(j)];
                cols=[cols;inLocs(i)];
            end
        end
        
        if ~isSymmetricDistFun
            d=distFun(x(inLocs(j):inLocs(j)+Lmin-1,:),x(inLocs(i):inLocs(i)+Lmin-1,:));
            if d<=th
                rows=[rows;inLocs(i)];
                cols=[cols;inLocs(j)];
            end
        end
    end
end
M=sparse(rows,cols,ones(numel(rows),1),T-Lmin+1,T-Lmin+1,numel(rows));
% clustering
clusters=clusterFun(M);
nC=numel(clusters);
toremove=[];
for c=1:nC
    if numel(clusters{c})<minSupport
        toremove=[toremove;c];
    else
        clusters{c}=sort(clusters{c});
        clusters{c}=clusters{c}(:);
    end
end
clusters(toremove)=[];
% convolution
step=1;
motifs=clusters;
maximalMotifs=cell(1,0);
for l=Lmin+1:step:Lmax
    cold=motifs;
    cnew=cell(1,0);
    nC=numel(cold);
    for i=1:nC
        ismaximal=isNotIncluded(cold{i},cnew,step);
        for j=1:nC
            f=directedIntersect(cold{i},cold{j},step);
            if ~isempty(f) && numel(f)>=minSupport
                [within,which]=included(f,cnew);
                if within
                    cnew{which}=chooseMaximal(f,cnew{which});
                else
                    cnew{numel(cnew)+1}=f;
                    %cnew=union(cnew,f);
                end
                if numel(f)==numel(cold{i})
                    ismaximal=false;
                end
            end
        end
        if ismaximal
            maximalMotifs{numel(maximalMotifs)+1}=[cold{i}(:),cold{i}(:)+l-2];
        end
    end
    if isempty(cnew)
        break;
    end
    if isempty(motifs)
        motifs=cell(1,0);
    end
    motifs=cnew;
end

locs =maximalMotifs(end:-1:1)';
if nargout>2
    stats=[];
    means=cell(1,0);
    nLoc=numel(locs);
    for i=1:nLoc
        [stats(i),means{i}]=calcStatsAndMean(locs{i},x);
    end
end

end

function maximal=chooseMaximal(c1,c2)
    maximal=union(c1,c2);
end
function notIncluded=isNotIncluded(short,long,step)
    notIncluded=true;
    if isempty(long)
        return;
    end
    n=numel(long);
    for i=1:n
        inter=intersect(short-step,long{i});
        if numel(inter)==numel(short)
            notIncluded=false;
            %break;
        end
        inter=intersect(short+step,long{i});
        if numel(inter)==numel(short)
            notIncluded=false;
            %break;
        end
    end
end
function d=distEuclidean(x1,x2)
    x1=x1(:); x2=x2(:);
    d=norm(x1-x2,2);
end
function clusters=findCliques(M)
    c=maximalCliques2(M);
    clusters=cell(size(c,2),1);
    for i=1:size(c,2)
        clusters{i}=find(c(:,i)~=0);
    end
end
function clusters=findConnectedComp(M)
    [labels]=graph_connected_components(M);
    n=max(labels);
    clusters=cell(n,1);
    for i=1:n
        clusters{i}=find(labels==i);
    end
end
function cplus=directedIntersect(c1,c2,step)
    cplus=[];
    nC1=numel(c1);
    for c=1:nC1
        if ~isempty(find(c2==c1(c)+step))
            cplus=[cplus;c1(c)];
        end
    end
end
function [within,which]=included(c,cplus) 
    within=false;
    which=0;
    if isempty(c) || isempty(cplus)
        return;
    end
    if iscell(cplus)
        n=numel(cplus);
        for i=1:n
            inter=intersect(c,cplus{i});
            if numel(inter)==numel(c) || numel(inter)==numel(cplus{i})
                within=true;
                which=i;
                break;
            end
        end
    else        
        inter=intersect(c,cplus);
        if numel(inter)==numel(c) || numel(inter)==numel(cplus)
            within=true;
        end
    end
end
function cunified=unify(c1,c2,cplus)
    within=false;
    cunified=[];
    if isempty(c) || isempty(cplus)
        return;
    end
    step=numel(cplus(1))-numel(c1(1));
    c12=directedIntersect(c1,c2,step);
    inter=intersect(c12,cplus);
    if numel(inter)==numel(c12) || numel(inter)==numel(cplus)
        within=true;
    end
    if within
        cunified=union(c12,cplus);
    end
end
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

function [locs,means,stats,centers,MDLs]=tanaka(x,th,nPointsPerSymbol,nSymbolsPerWord,...
    alphabetSize,Lmin,Lmax,minSupport,distFun,removeSimilarMotifs,maxOverlap)
% Implements Tanaka et al.'s multidimensional time-series motif discovery method
[T,n]=size(x);
if ~exist('maxOverlap','var') || isempty(maxOverlap)
    maxOverlap=0.1;
end
if ~exist('removeSimilarMotifs','var') || isempty(removeSimilarMotifs)
    removeSimilarMotifs=1;
end
if ~exist('Lmax','var') || isempty(Lmax)
    Lmax=T-Lmin;
end
if ~exist('distFun','var') || isempty(distFun)
    distFun=@dEN;
end
if ~exist('minSupport','var') || isempty(minSupport)
    minSupport=2;
end
%nBSPerMotif=Lmin/(nPointsPerSymbol*nSymbolsPerWord);

if n>1
x=tspca(x);
end
% % 
[xd]=timeseries2symbol(x,nPointsPerSymbol*nSymbolsPerWord,...
    nSymbolsPerWord,alphabetSize,1);
[symb,s2bs]=generateBS(xd);
sa=size(s2bs,1);


MDLs=[];
locs=cell(1,0);
%Ls=[];
centers=cell(1,0);
for L=Lmin:Lmax
    nBSPerMotif=L-(nPointsPerSymbol*nSymbolsPerWord)+1;
    md=ts2mat(symb,nBSPerMotif);
    [bss,s2y]=generateBS(md);
    [counts,loc]=countN(s2y,bss);
    [scounts,indx]=sort(counts,'descend');
    sloc=loc(indx);

    nBSs=numel(counts);
    for i=1:nBSs
        bs=s2y(bss(indx(i)),:);
        bsLoc=sloc{i};
        c=scounts(i);
        if c<minSupport
            break;
        end
        A=zeros(c,c);
        tsLoc=bsLoc;
        for j=1:c
            for k=i+1:c
                A(j,k)=distFun(x(tsLoc(j):tsLoc(j)+L-1),x(tsLoc(k):tsLoc(k)+L-1));
                A(k,j)=A(j,k);
            end
        end
        B=A<=th;
        for j=1:c
            B(j,j)=0;
        end
        bcount=sum(B,2);
        [sbcount,bidx]=sort(bcount,'descend');
        sbcount(sbcount<max(sbcount))=0;
        occurrences=bidx(sbcount>0.5);
        dists=sum(A(occurrences,:),2);

        [dcount,didx]=sort(dists,'ascend');
        centers{numel(centers)+1}=bidx(didx(1));
        %Ls(numel(Ls)+1)=L;
        
        locs{numel(locs)+1}=[sloc{didx(1)},sloc{didx(1)}+L-1];
        MDLs(numel(MDLs)+1)=mdl2(xd,bs,sa,numel(locs{end}));
    end
end

[MDLs,indx]=sort(MDLs,'ascend');
centers=centers(indx)';
locs=locs(indx)';

[locs,toremove]=mdRemoveNonMaximal(locs);
MDLs(toremove)=[];
centers(toremove)=[];


if removeSimilarMotifs
    [locs,centers,MDLs]=removeNearMotifs(x,locs,centers,distFun,th, ...
        MDLs,sa,nPointsPerSymbol,nSymbolsPerWord,maxOverlap);
    [MDLs,indx]=sort(MDLs,'ascend');
    centers=centers(indx)';
    locs=locs(indx)';
end

if nargout>2
    stats=[];
    means=cell(1,0);
    nLoc=numel(locs);
    for i=1:nLoc
        [stats(i),means{i}]=calcStatsAndMean(locs{i},x);
    end
end

end

function [locs,centers,MDLs,toremove]=removeNearMotifs(x,locs,centers,distFun...
    ,tau,MDLs,sa,nPointsPerSymbol,nSymbolsPerWord,maxOverlap)
[T,n]=size(x);
C=numel(centers);
toremove=[];
for i=1:C
    if isempty(locs{i})
        continue;
    end
    L1=locs{i}(1,2)-locs{i}(1,1)+1;
    for j=i+1:C
        if isempty(locs{j})
            continue;
        end
        L2=locs{j}(1,2)-locs{j}(1,1)+1;
        if L2~=L1
            continue;
        end
        if distFun(x(centers{i}:centers{i}+L1-1),x(centers{j}:centers{j}+L2-1))<2*tau
            K=size(locs{j},1);
            for k=1:K
                if distFun(x(centers{i}:centers{i}+L1-1),x(locs{j}(k,1):locs{j}(k,2)))<tau
                    Ki=size(locs{i},1);
                    overlap=0;
                    for l=1:Ki
                        overlap=max([overlap,overlapfraction(locs{i}(l,:),locs{j}(k,:))]);
                    end
                    if overlap<maxOverlap
                        locs{i}(end+1,:)=locs{j}(k,:);
                        MDLs(i)=updatemdl2(T-L1+1,L1-nPointsPerSymbol*nSymbolsPerWord+1,sa,size(locs{i},1)-1,MDLs(i));
                    end
                end
            end
            toremove=[toremove j];
        end
    end
end
locs(toremove)=[];
centers(toremove)=[];
MDLs(toremove)=[];
end

function locs =findLocs(bss,bs)
    locs=find(sum(abs(bss-repmat(bs,T,1)),2)<0.01);
end
function d=dEN(x1,x2)
    x1=x1-mean(x1);
    x2=x2-mean(x2);
    d=norm(x1-x2,2)./numel(x1);
end

function [counts,locs]=countN(s2x,xd)
    %counts=hist(s2x,size(xd,1));
    ns=size(s2x,1);
    counts=zeros(ns,1);
    locs=cell(ns,1);
    for i=1:ns
        locs{i}=find(xd==i);
        counts(i)=numel(locs{i});
    end
end

function m=ts2mat(x,l)
T=numel(x);
m=zeros(T-l+1,l);
for i=1:T-l+1
    m(i,:)=x(i:i+l-1);
end
end

function [y,s2x]=generateBS(x)
    [s2x,ia,y]=unique(x,'rows','stable');
    %s2x=s(ia,:);
    
end


function v=mdl(x,sc,sa)
T=numel(x);
np=numel(sc);
sp=numel(unique(sc));

q=nOccurrences(x,sc);
na=T-(np-1)*q;

v=log2(np)+np*log2(sp)+log2(na)+na*log2(sa+q);

end

function v=mdl2(x,sc,sa,q)
T=numel(x);
np=numel(sc);
sp=numel(unique(sc));

na=T-(np-1)*q;

v=log2(np)+np*log2(sp)+log2(na)+na*log2(sa+q);

end

function v=updatemdl2(T,np,sa,q,oldv)
na=T-(np-1)*q;
nanew=T-(np-1)*(q+1);
v=oldv-log2(na)+log2(nanew)...
    -na*log2(sa+q)+nanew*log2(sa+q+1);

end

function n=nOccurrences(x,sc)
T=numel(x);
m=numel(sc);
n=0;
for t=1:T-m+1
    if sum(abs(x(t:t+m-1)-sc))<0.0001
        n=n+1;
    end
end
end
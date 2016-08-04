function [begends,isChanged] = detectMotifDefault(x,begends,meanMotif,...
    stats,candLocs,maxSmallDist,minLargeDist,dFun,dFunParams...
    ,stopFun,stopFunParams,optional)
%Detects all occurrences of a motif in a time series
%
% Inputs:
% =======
% begends       a n*2 array giving the beginning and end of already found
%               motif occurrences
% meanMotif     The mean of the found occurrences
% stats         statistics of the found occurrences (see DGR.m for its
%               contents)
% candLocs      Locations at which it is likely to find motifs. It may be
%               empty
% optional      a cell array that is used to pass empty parameters. The
%               function must work if this is passed as empty or not even
%               passed at all
% 
% Outputs:
% =======
% newbegends    The locations of all found motif occurrences  in the same 
%               format as begends
% isChanged     nonzero if newbegends is not the same as begends
%
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%

T=length(x);
x=x(:)';
stepSize=1;

len=length(meanMotif);
if ~isempty(candLocs)
    nL=length(candLocs);
    locs=[];
    for i=1:nL
        locs=[locs,(-len:stepSize:+len)+candLocs(i)];
    end
    locs(locs<1)=[];
    locs(locs>T-len)=[];
    locs=unique(locs);  
else
    locs=1:stepSize:T-len;
end

nL=length(locs);
toreuse=[];
savedbe=begends;
toignore=zeros(nL,1);
nB=size(begends,1);
for i=1:nB
    l=(begends(i,2)-begends(i,1)+1);
    for j=1:nL        
        if locs(j)>=begends(i,1)-l+1 && locs(j)<=begends(i,2)
            toignore(j)=1;
        end
    end
end
for i=1:nL
    if toignore(i); continue; end;
    d=dFun(meanMotif,x(locs(i):locs(i)+len-1),dFunParams);    
    if d>maxSmallDist; continue; end;
    newbegends=[begends;locs(i),locs(i)+len-1];
    if size(begends,1)>1
        [stop,toreuse]=stopFun(x,begends,newbegends,toreuse,stopFunParams);
        if stop; continue; end;
    end
    begends=newbegends;
    toignore(locs<newbegends(end,2))=1;    
end
isChanged=(size(savedbe,1)<size(begends,1));   
end


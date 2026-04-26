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

function [locs,means,stats] = mcfull(x,lrange,candLocs,varargin)
%Modified catlano algorithm. implements mcfull algorithm and if no candLocs
%are given 
%
%
% Inputs:
% =======
% x         Input time series (can be 2-dimensional but the distance
%           function should be able to find distances between the multiple dimensions
%           together. T (the length of x) is size(x,2) which means that for
%           single dimensional data x must be a ROW vector. if x is
%           multidimensional (m*T) then tspca will be used to convert it
%           into a single dimensional array before doing anything else.
% lrange    The minimum and maximum expected lengths of motifs. If a
%           single number is passed it is treated as the maximim and the
%           minimum is automatically calcualted from it as:
%           (max(3,floor(0.9*lrange))
% candLocs  [OPTIONAL] Candidate locations of motifs found from the constraints. The
%           simplest method to get these is to RSST the signal and localize
%           the changes.
% varargin  See the switch statement for their meaning
%
% Outputs:
% ========
% locs      a m*1 cell array representing the locations of m motifs in the
%           time series. each element of the cell-array contains a n_i*2
%           array representing the beginning and end of one occurrence of
%           the motif
% means     The means of the motifs found as a m*1 cell array each with a
%           mean of one motif
% stats     the statistics of each motif in the form of a structure array.
%           members:
%               mean    The mean difference between corresponding
%                           points in the occurrences
%               var     The standard deviation difference between corresponding
%                           points in the occurrences
%               max     The maximum difference between corresponding
%                           points icdn the occurrences
%               min     The minimum difference between corresponding
%                           points in the occurrences
%               median  The mmedian difference between corresponding
%                           points in the occurrences
%               mode    The mode of the difference between corresponding
%                           points in the occurrences
%
% Notes
% =====
% - either candLocs or nWindows must be given. if both
%   are given, chgScores will be ignored
%
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, Constrained Motif Discovery in Time
% Series, New Generation Computing, 27(2009)319-346 
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%
%
T=length(x);
if length(lrange)<2
    lrange=[max(3,floor(0.9*lrange)),lrange];
end
if nargin<3
    candLocs=[];    
end

alpha=1;   % based in results in orginial catalano's paper
K=[];       % will be set to nWindows-1 (the maximum possible later)
chgScores=[];
nWindows=10; % based in results in orginial catalano's paper


wbar=[];
nSkip=[];
extraLength=0.2;
nCompWindows=[];  % default is 0.5nWindows
dFun=@dE;
dFunParams=[];
stopFun=@stopIfAboveMean;
stopFunParams=[];
dPointDiffFun=@abs;     % will be used only if we need to return stats
doDetect=false;

nArgs=size(varargin,2);
if(nArgs==1 && isempty(varargin{1}))
   nArgs=0;
elseif(nArgs==1 && iscell(varargin{1}))
    varargin=varargin{1};
    nArgs=size(varargin,2);
end
if(nArgs>0)
    if(mod(nArgs,2)~=0)
        error('The optional arguments must be in the form name,value so they must be even!!!');
    end
    for i=1:2:nArgs
        switch(lower(varargin{i}))                        
            case {'nwindows','nw'} % [OPTIONAL] The number of candidate windows to use. Default is 10                                    
                nWindows=(varargin{i+1});
            case {'chgscores','chg'}     % [OPTIONAL] a time series of the same length as x. the score
                                         %             at each point is proportional to finding a motif occurrence 
                                         %             around this point. The easiest way to get this score is
                                         %             through a change point discovery routine like rsst() or cpd()
                chgScores=(varargin{i+1});
            case {'k','neffectivedists'} %[OPTIONAL] the number of candidate windows nearest to each
                                         %           subsequence to keep for further processing. Default =nWindows-1
                K=(varargin{i+1});
            case {'alpha','distdiscount'}%[OPTIONAL]significance used for removing noise-like candidate
                                         %           motifs
                alpha=(varargin{i+1});
            otherwise
                error('unknown command: %s',varargin{i});
        end
    end
end

% initialize locations around which we will sample
if ~isempty(candLocs)
    nWindows=length(candLocs);
else
    if isempty(chgScores)
        if nWindows==0
            error('you cannot give neither nWinodws nor candLocs');
        end
        candLocs=randi(T,nWindows,1);
    else
        if nWindows==0
            error('you cannot give chgScores without giving nWindows>0');
        end
        candLocs=zeros(nWindows,1);
        %normalize chgScores to be a probability distribution
        chgScores=chgScores./((chgScores(:))'*(chgScores(:)));
        cumchgscore=cumsum(chgScores);
        for i=1:nWindows
            p=rand(1,1);
            for j=1:T
                if p>=cumchgscore(j)
                    caldLocs(i)=j;
                    break;
                end
            end
        end
    end
end

% set parameters
w=round((1+extraLength)*lrange(2));
if isempty(wbar)
    wbar=max(3,ceil(0.1*w));
end
if isempty(nSkip)
    nSkip=max(1,floor(.5*wbar));
end

nSub=numel([0:nSkip:w-wbar]);

% generate the noise window
noise=x(randi(T,w,1));

%initialize other variables
if isempty(nCompWindows)
    nCompWindows=max(2,ceil(0.5*nWindows));
end
iterations=ceil(nCompWindows/nWindows);
if isempty(K)
    K=nWindows-1;
else
    K=min(K,nWindows-1);
end

[allSW]=createSWSet(candLocs,nWindows,nSub,w,T,nSkip);
[noiseSW]=createSWSet(ceil(w/2),1,nSub,w,w,nSkip);

locs=[];
means=[];
stats=[];

%find nearest subsequences to noise sequences    
[bestMatchesMn]=compareFindBests...
    (noise,x,allSW,noiseSW,wbar,K,dFun,dFunParams);
maxSmallDist=mean(bestMatchesMn(:,1));
gamma=alpha*maxSmallDist;
% find nearest subsequences to candidate subsequences

activeSubs=ones(nSub,1);
for currentWin=1:nWindows        
    for iteration=1:iterations
        if ~any(activeSubs)
            break;
        end
        candSW=allSW(currentWin,activeSubs~=0);        
        [bestMatchesM]=compareFindBests...
            (x,x,allSW(randi(nWindows,nCompWindows,1),:),...
            candSW,wbar,K,dFun,dFunParams);
        [activeSubs]=updateCandSW(activeSubs,bestMatchesM,gamma);
    end
    if ~any(activeSubs)
        continue;
    end
    mySubs=allSW(currentWin,:);
    begend=[-1,-1];
    nNotMatched=0;
    current=1;
    for i=1:nSub
        if begend(current,1)<0
            if activeSubs(i)
                begend(current,1)=allSW(currentWin,i);
                begend(current,2)=allSW(currentWin,i)+wbar-1;                
            end
        else
            if activeSubs(i)
                begend(current,2)=allSW(currentWin,i)+wbar-1;
                nNotMatched=0;
            else
                nNotMatched=nNotMatched+nSkip;
                if nNotMatched>=wbar
                    begend=[begend;-1,-1];
                    current=current+1;
                    nNotMatched=0;
                end
            end
        end
    end
    begend(begend(:,1)<0,:)=[];
    if doDetect   
        [dummy,iMax]=max(begend(:,2)-begend(:,1));
        begend=begend(iMax,:);
        means=[means,{x(begend(1):begend(2))}];
    
        begend=detectMotifDefault(x,begend,means{end},...
            [],candLocs,gamma,maxSmallDist,dFun,dFunParams...
            ,stopFun,stopFunParams);
    else
        means=[means,{x(begend(1,1):begend(1,2))}];
    end
    if size(begend,1)>1
        locs=[locs;{begend}];
    else
        means(end)=[];
    end
end

if nargout>2
    nMotifs=length(locs);
    for i=1:nMotifs
        [stats,means]=updateStatsAndMean(locs{i},x,dPointDiffFun,stats,means,i);
    end
end
end
function [activeSubs]=updateCandSW(activeSubs,bestMatchesM,gamma)
        %,bestMatchesD,bestMatchesW,bestMatchesS)    
    nSub=length(activeSubs);    
    j=1;
    for k=1:nSub
        if ~activeSubs(k)
            continue;
        end
        if bestMatchesM(j)>gamma
            activeSubs(k)=0;
            j=j+1;
            continue;
        end
        %matchSet{k}=[matchSet{k};bestMatchesW(j,:),bestMatchesS(j,:)];
        %matchSetD{k}=[matchSet{k};bestMatchesD(j,:)];
        j=j+1;
    end    
end
function [bestMatchesM,bestMatchesD,bestMatchesW,bestMatchesS]=compareFindBests(...
    cand,comp,allSW,candSW,wbar,K,dFun,dFunParams)    
    nWindows=size(allSW,1);
    nSub=size(allSW,2);
    nSubSrc=length(candSW);
    dists=zeros(nWindows,nSub);    
    dists(:,:)=inf;    
    bestMatchesW=zeros(nSubSrc,K);
    bestMatchesS=zeros(nSubSrc,K);
    bestMatchesD=zeros(nSubSrc,K);
    bestMatchesM=zeros(nSubSrc,1);
        
    for k=1:nSubSrc
        matches=zeros(nWindows);
        matchDists=zeros(nWindows);        
        for j=1:nWindows            
            for l=1:nSub
                dists(j,l)=dFun(comp(allSW(j,l):allSW(j,l)+wbar-1),...
                    cand(candSW(1,k):candSW(1,k)+wbar-1),dFunParams);                
            end           
            [matchDists(j),matches(j)]=min(squeeze(dists(j,:)));
        end
        [bestMatchesAllD,bestMatchesAllW]=sort(matchDists,'ascend');
        bestMatchesD(k,:)=bestMatchesAllD(1:K);
        bestMatchesW(k,:)=bestMatchesAllW(1:K);
        bestMatchesS(k,:)=matches(bestMatchesW(k,:));        
        bestMatchesM(k)=mean(bestMatchesD(k,:));
    end    
    [bestMatchesM,I]=sort(bestMatchesM,'ascend');
    if nargout>1
        bestMatchesD=bestMatchesD(I,:);
        bestMatchesS=bestMatchesS(I,:);
        bestMatchesW=bestMatchesW(I,:);
    end
end

function [loc]=createSWSet(candLoc,nWindows,nSub,w,T,nSkip)
% creates and returns a subwindow set in the form of a nWindows*nSub*2 matrix
loc=zeros(nWindows,nSub);
for i=1:nWindows
    begend=min(T-w+1,max(1,candLoc(i)-floor(w/2)));        
    for j=1:nSub
        loc(i,j)=begend+(j-1).*nSkip;        
    end
end
end

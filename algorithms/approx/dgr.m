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

function [locs,means,stats,maxSmallDist,minLargeDist] = dgr(x,lrange,candLocs,varargin)
%Distance Graph Relaxation CMD algorithm
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
% candLocs  Candidate locations of motifs found from the constraints. The
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
%                           points in the occurrences
%               min     The minimum difference between corresponding
%                           points in the occurrences
%               median  The mmedian difference between corresponding
%                           points in the occurrences
%               mode    The mode of the difference between corresponding
%                           points in the occurrences
%
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad, Toyoaki Nishida, and Shogo Okada, Unsupervised
% Simultaneous Learning of Gestures, Actions and their Associations for 
% Human-Robot Interaction, IEEE IROS 2009 (October 11 to 15, St. Louis, 
% MO, USA), pp. 2537-2544
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%

% process inputs
if length(lrange)<2
    lrange=[max(3,floor(0.9*lrange)),lrange];
end

% intialize internal parameters
stemLen=[];
maxLen=[];
minLen=[];
aroundRange=[];
expectedNMotifs=[];  
maxStemsToTry=[];
maxStemLocs=[];
maxStemsAlreadyUsedConsecutive=inf;
epsilon=[];
acceptableMeanIncrease=[];
nRandCandLocs=0;

visualize=0;
outputFileID=1;
verbos=0;
minExtension=[];             % minimum number of points to extend 

tryExtendingStemsSlowly=0;
useBisection=1;
tryCombiningNearStems=1;
tryCombiningSingeltonUnusedStems=0;
tryGrowingSingeltonUnusedStems=1;

useMinMaxDists=1;
significance=0.05;
ignoreOverlapsDuringExtension=1;
ignoreOverlapsDuringCombniation=1;
epsilonDistance=1e-8;

dPointDiffFun=@abs;
dFun=@dE;
dFunParams=[];
stopFun=@stopIfAboveMean;
stopFunParams=[];
optimalDistFun=@optimalDistKMeans; %[]
optimalDistFunParams=[];
% change internal parameters using optional arguments

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
            case {'trycombiningnearstems','cns'} % whether or not stems that look similar ar to be combined in a single motif stem during first stage
                tryCombiningNearStems=(varargin{i+1});   
            case {'tryextendingstemsslowly','ess'} % whether or not the stems are extended slowly by adding small stems to it
                tryExtendingStemsSlowly=(varargin{i+1});               
            case {'trycombiningsingeltonunusedstems'}          
                tryCombiningSingeltonUnusedStems=(varargin{i+1});                   
            case {'trygrowingsingeltonunusedstems'}          
                tryGrowingSingeltonUnusedStems=(varargin{i+1});                   
            case {'usebisection','bisect'}          % whether or not to use bisection when trying to fastly extend stems
                useBisection=(varargin{i+1});                   
            case {'ignoreoverlapsduringextension'}          
                ignoreOverlapsDuringExtension=(varargin{i+1});                   
            case {'ignoreoverlapsduringcombniation'}          
                ignoreOverlapsDuringCombniation=(varargin{i+1});                   
            case {'dpointdifffun'}          
                dPointDiffFun=(varargin{i+1});                   
            case {'dfun'}          
                dFun=(varargin{i+1});                   
            case {'dfunparams'}          
                dFunParams=(varargin{i+1});                   
            case {'stopfun'}          
                stopFun=(varargin{i+1});                   
            case {'stopfunparams'}          
                stopFunParams=(varargin{i+1});                   
            case {'optimaldistfun'}          
                optimalDistFun=(varargin{i+1});                   
            case {'optimaldistfunparams'}          
                optimalDistFunParams=(varargin{i+1});                       
            case {'acceptablemeanincrease'}          
                acceptableMeanIncrease=(varargin{i+1});                   
            case {'nrandcandlocs'}          
                nRandCandLocs=(varargin{i+1});
            case {'visualize'}          
                visualize=(varargin{i+1});                   
            case {'outputfileid'}          
                outputFileID=(varargin{i+1});                   
            case {'verbos'}          
                verbos=(varargin{i+1});                   
            case {'minextension'}          
                minExtension=(varargin{i+1});              
            case {'useminmaxdists'}          
                useMinMaxDists=(varargin{i+1});                   
            case {'significance'}          
                significance=(varargin{i+1});                               
            case {'epsilondistance'}          
                epsilonDistance=(varargin{i+1});                                               
            case {'stemlen','sl'}          
                stemLen=(varargin{i+1});               
            case {'expectednmotifs'}          
                expectedNMotifs=(varargin{i+1});               
            case {'maxstemstotry'}          
                maxStemsToTry=(varargin{i+1});                   
            case {'maxstemlocs'}          
                maxStemLocs=(varargin{i+1});                   
            case {'maxstemsalreadyusedconsecutive'}          
                maxStemsAlreadyUsedConsecutive=(varargin{i+1});                   
            case {'epsilon'}          
                epsilon=(varargin{i+1});                               
            case {'maxlen'}          
                maxLen=(varargin{i+1});                       
            case {'minlen'}          
                minLen=(varargin{i+1});                               
            case {'aroundrange'}          
                aroundRange=(varargin{i+1});                       
            otherwise
                error('unknown command: %s',varargin{i});
        end
    end
end

% recalculate any needed parameters;
if verbos
    begninning=tic;
end
if ~isvector(x)
    if size(x,1)>size(x,2)
        x=x';
    end
    x=tspca(x);
else
    if size(x,1)~=1
        x=x';
    end
end
T=length(x);
if isempty(epsilon)
    epsilon=0.000001*(max(x)-min(x));
end
if isempty(acceptableMeanIncrease)
    acceptableMeanIncrease=2*significance;
end
if isempty(maxStemsToTry)
    maxStemsToTry=T;
end
if isempty(stemLen)
    stemLen=floor(lrange(1)/2);
end
if isempty(maxLen)
    maxLen=lrange(2);
end
if isempty(minLen)
    minLen=lrange(1);
end
if isempty(nRandCandLocs)
    nRandCandLocs=ceil(log2(T));
end
if isempty(aroundRange)
    aroundRange=[0,-stemLen,round([-maxLen:stemLen:maxLen])];
end
nCandidatesPerLoc=length(aroundRange);
if isempty(stopFunParams)
    stopFunParams={significance,epsilon,dPointDiffFun,acceptableMeanIncrease};
end
if isempty(expectedNMotifs)
    expectedNMotifs=length(candLocs);  
end
if isempty(maxStemLocs)
    maxStemLocs=max(3,round(T*minLen/(mean(lrange)*stemLen)));
end
if isempty(minExtension)
    minExtension=1;
end
% % remove very near change locations
% candLocs=sort(candLocs);
% candLocDist=candLocs(2:end)-candLocs(1:end-1);
% near=find(candLocDist<2*lrange(2));
% while(~isempty(near))
%     candLocs(near(max(1,floor(length(near)/2))))=[];
%     candLocDist=candLocs(2:end)-candLocs(1:end-1);
%     near=find(candLocDist<2*lrange(2));
% end

% inititalize
if nRandCandLocs>0
    candLocs=[candLocs(:);randi(T-stemLen,nRandCandLocs,1)];
end

locs=[];
means=[];
savedStats=[];
nLocs=length(candLocs);


% find stem locations around candLocs
if verbos; fprintf(outputFileID,'calculating candidate locations ... '); start=tic; end;
stemLocs=zeros(nLocs,nCandidatesPerLoc);   % this will give the BEGINNING of the stems we will use in the distance graph
for loc=1:nLocs
    tmp=aroundRange+candLocs(loc);
    if(tmp(1)<1); 
        tmp=tmp-tmp(1)+1; 
    end;
    if(tmp(end)>T-stemLen+1); 
        tmp=tmp-(tmp(end)-T+stemLen-1); 
    end;    
    stemLocs(loc,:)=tmp;
end
stemLocs=reshape(stemLocs,numel(stemLocs),1);
stemLocs(stemLocs<1)=[];
stemLocs(stemLocs>T)=[];
stemLocs=unique(stemLocs);
if length(stemLocs)>3
    toremove=zeros(size(stemLocs));
    for i=1:length(stemLocs)-1    
        if toremove(i); continue; end;
        for j=i+1:length(stemLocs)
            if(stemLocs(j)-stemLocs(i)<minLen)
                toremove(j)=1;
            else            
                break;
            end
        end
    end
    stemLocs(toremove~=0)=[];
end
stemLocs=sort(stemLocs);
nstems=length(stemLocs);
if(nstems>maxStemLocs)
    toremove=round([1:nstems/maxStemLocs:nstems]);
    stemLocs=stemLocs(toremove);
end
nstems=length(stemLocs);
if nstems<3
    warning('DGR::NoEnoughStems',sprintf('only %d stems. we need at least 3',nstems));
    locs=[];
    means=[];
    stats=[];
    return;
end
ndists=nstems*(nstems-1)/2;
if verbos; fprintf(outputFileID,'Done in %gs [nstems=%d,stemLen=%d,percentOfT=%05.3f]\n'...
        ,toc(start),nstems,stemLen,double(nstems*stemLen)/double(T)); end;
% now we have the stem locations and we build the distance matrix
if verbos; fprintf(outputFileID,'finding all distances and arranging them ... '); start=tic; end;
% distStats=zeros(nstems,3);
if tryCombiningNearStems
    dists=zeros(nstems,nstems);
end
distLinear=zeros(ndists,6);
nxt=1;
for i=1:nstems-1
    for j=i+1:nstems
        dd=dFun(x(stemLocs(i):stemLocs(i)+stemLen-1),...
                        x(stemLocs(j):stemLocs(j)+stemLen-1),dFunParams);
        if tryCombiningNearStems
            dists(i,j)=dd;
        end
        distLinear(nxt,:)=[stemLocs(j),stemLocs(i),round(dd./epsilonDistance)...
            ,j,i,(stemLocs(j)-stemLocs(i))+1];
        nxt=nxt+1;
    end
%     distStats(i,1)=median(dists(i,i+1:end));
%     distStats(i,2)=mean(dists(i,i+1:end));
%     distStats(i,3)=std(dists(i,i+1:end));
end
if tryCombiningNearStems
    dists=dists./max(max(dists));
end

%[dummy,indx]=sort(distLinear(:,3),'ascend');
%distLinear=distLinear(indx,:);
distLinear=sortrows(distLinear,[3,-6]);
distLinear(:,3)=distLinear(:,3)*epsilonDistance;
if verbos; fprintf(outputFileID,'Done in %gs\n',toc(start)); end;
%calculate the number of stems to use
if verbos; fprintf(outputFileID,'Finding optimal distance ... '); start=tic; end;
if ~isempty(optimalDistFun)
    [nStemsToTry,maxSmallDist,minLargeDist]=optimalDistFun(distLinear,optimalDistFunParams);
else
    nStemsToTry=expectedNMotifs*2;
    maxSmallDist=distLinear(nStemsToTry,3);
    minLargeDist=distLinear(min(2*nStemsToTry,size(distLinear,1)*.8),3);
end
if ~useMinMaxDists
    maxSmallDist=inf; minLargeDist=inf;
end
nStemsToTry=min(nStemsToTry,maxStemsToTry);
if verbos; fprintf(outputFileID,'Done in %gs [nStemsToTry=%d of %d,maxSmall=%g,minLarge=%g]\n',...
        nStemsToTry,size(distLinear,1),maxSmallDist,minLargeDist,toc(start)); end;
% [distdiff,locMaxDistDif]=sort(thin((distLinear(2:end,3)-distLinear(1:end-1,3))./(distLinear(2:end,3)+distLinear(1:end-1,3))),'descend');
% locMaxDistDif(isinf(distdiff))=[]; locMaxDistDif(isnan(distdiff))=[];
% distdiff(isinf(distdiff))=[]; distdiff(isnan(distdiff))=[];
% distLimits=distLinear(locMaxDistDif);
% if(distLimits(1)<1e-20); distLimits(1)=[]; end;


% for all distances found starting from smallest to largest
stemUsed=zeros(size(stemLocs));
stemTryed=[];
if verbos; fprintf(outputFileID,'Testing stems:\n'); ss=tic; end;
if isempty(maxStemsAlreadyUsedConsecutive)
    maxStemsAlreadyUsedConsecutive=0.125*nStemsToTry;
end
stemsBothUsed=0;
if visualize
    f=figure;    
end
for i=1:nStemsToTry     
    if verbos>1; fprintf(outputFileID,'\t[%d/%d] Combination Trial ... ',i,nStemsToTry); start=tic; end;    
    stemTryedAround=zeros(size(stemLocs));
    stemTryed{i}=zeros(size(stemLocs));
    if(stemUsed(distLinear(i,4)) && stemUsed(distLinear(i,5)))        
        if verbos>1; fprintf(outputFileID,'Done in %gs[continued1]\n',toc(start));   end;    
        stemsBothUsed=stemsBothUsed+1;
        if stemsBothUsed>= maxStemsAlreadyUsedConsecutive
            break;
        else
            continue;
        end        
    elseif (stemUsed(distLinear(i,4)) || stemUsed(distLinear(i,5)))
        stemsBothUsed=0;
        if ~tryCombiningSingeltonUnusedStems
            if verbos>1; fprintf(outputFileID,'Done in %gs[continued2]\n',toc(start));   end;    
            if ~tryGrowingSingeltonUnusedStems
                continue;
            end
        elseif ~isempty(locs)
            unused=5-stemUsed(distLinear(i,5));
            for m=1:length(locs)
                stemTryed{m}(distLinear(i,unused))=1;
                [stop,isoverlapping,stemUsed,locs{m},means{m}]=tryCombining...
                    (distLinear(i,unused),stemLocs,...
                    stemUsed,locs{m},means{m},savedStats{m},length(means{m}),x,T,stemLen...
                    ,dFun,dFunParams,stopFun,stopFunParams,0,ignoreOverlapsDuringCombniation);
                
            end            
            if verbos>1; fprintf(outputFileID,'Done in %gs\n',toc(start));   end;    
            continue;
        else
            if verbos>1; fprintf(outputFileID,'Done in %gs[continued3]\n',toc(start));   end;    
            continue;
        end
    else
        stemsBothUsed=0;
    end
    if verbos>1; fprintf(outputFileID,'Done in %gs\n',toc(start));   end;    
    if verbos>1; fprintf(outputFileID,'\t[%d/%d] Main ... :\n',i,nStemsToTry); start=tic; end;        
    begends=[distLinear(i,1:2)',distLinear(i,1:2)'+stemLen-1];
    stemUsed(distLinear(i,4:5))=1;
    % slid the second stem member until it gives minimum distance to the
    % first    
    if( size(begends,1)==2)
        mind=distLinear(i,3); minj=0;
        for j=[-stemLen:-1,1:stemLen]
            if(begends(2,1)+j>0 && begends(2,2)+j<T+1)
                dtmp=dFun(x(begends(1,1):begends(1,2)),...
                                x([begends(2,1):begends(2,2)]+j),dFunParams);
                if(dtmp<mind)
                  mind=dtmp; minj=j;
                end
            end
        end
        begends(2,:)=begends(2,:)+minj;
    end
    toreuse=[];
    if(visualize); showbegends(x,T,begends,1,maxLen); end;        
    
    if verbos>2; fprintf(outputFileID,'\t\t[%d of %d] Extension ...',i,nStemsToTry); s1=tic; end;    
    % try to extend the stem one motif per step to both sides then left
    % then right then by a single point    
    if tryExtendingStemsSlowly
        [begends,toreuse]=extendStem([-stemLen,stemLen],x,T,maxLen,begends ...
            ,stopFun,[],stopFunParams,maxSmallDist,minLargeDist,dFun,...
            dFunParams,minExtension,tryExtendingStemsSlowly,useBisection...
            ,0,ignoreOverlapsDuringExtension);
        if(visualize); showbegends(x,T,begends,0,0); end;        
    else
        [begends,toreuse]=extendStem([-maxLen+stemLen,0],x,T,maxLen,begends ...
            ,stopFun,[],stopFunParams,maxSmallDist,minLargeDist,dFun,dFunParams...
            ,minExtension,tryExtendingStemsSlowly,useBisection,0,ignoreOverlapsDuringExtension);
        if(visualize); showbegends(x,T,begends,0,0); end;
        [begends,toreuse]=extendStem([0,maxLen-stemLen],x,T,maxLen,begends ...
            ,stopFun,toreuse,stopFunParams,maxSmallDist,minLargeDist,dFun...
            ,dFunParams,minExtension,tryExtendingStemsSlowly,useBisection...
            ,0,ignoreOverlapsDuringExtension);
        if(visualize); showbegends(x,T,begends,0,0); end;
    end
    if verbos>2; fprintf(outputFileID,' Done %gs\n',toc(s1)); end;
    if verbos>2; fprintf(outputFileID,'\t\t[%d of %d] Combination ...',i,nStemsToTry); s2=tic; end;
    % find the mean of this motif
    mLen=begends(:,2)-begends(:,1)+1;
    minoLen=min(mLen);  %maxoLen=max(mLen);
    meanMotif=zeros(1,minoLen);
    for kk=1:2
        meanMotif=meanMotif+(x(begends(kk,1):begends(kk,2)));
    end
    meanMotif=meanMotif./2;
    % check if any one of the stems we still did not test is within the
    % new extended range
    tryCombine=tryCombiningNearStems;
    while(tryCombine)
        tryCombine=0;                
        for s=1:nstems            
            if(stemTryedAround(s))
                continue;
            end
            nStemsInGroup=size(begends,2);
            for k=1:nStemsInGroup
                if stemLocs(s)>=begends(k,1) && stemLocs(s)<=begends(k,2)
                    stemTryedAround(s)=1;
                    % try to combine motifs nearest to s in the motif
                    d=[[s+1:nstems,1:s-1]',[dists(s,s+1:end),dists(1:s-1,s)']' ];
                    [dummy,index]=sort(d(:,2),'ascend');
                    d=d(index,:);
                    maxNear=find(d(:,2)-maxSmallDist,1,'first');
                    nTrials=min(nStemsToTry,maxNear);
                    for near=1:nTrials;
                        % confirm that this one is not included in
                        % something else
                        current=d(near,1);
                        if(stemTryed{i}(current) || stemUsed(current))
                            continue;
                        end
                        stemTryed{i}(current)=1;
                        [stop,isoverlapping,stemUsed,begends,meanMotif]= ...
                            tryCombining(current,stemLocs,...
                            stemUsed,begends,meanMotif,toreuse,minoLen,x,T...
                            ,stemLen,dFun,dFunParams,stopFun,stopFunParams,0 ...
                            ,ignoreOverlapsDuringCombniation);
                        if isoverlapping; 
                            stemUsed(current)=1;
                            continue; 
                        end;
                        if stop                             
                            continue;
                        else
                            tryCombine=1;
                            if(visualize); showbegends(x,T,begends,0,0); end; 
                        end;
                    end
                    break;
                end
            end                             
        end
    end
    dst=begends(:,2)-begends(:,1)+1;
    if max(dst)>=minLen
        locs=[locs;{begends}];    
        means=[means;{meanMotif}];        
        savedStats=[savedStats;{toreuse}];
    end
%     if visualize
%         close(f);
%     end
    if verbos>1; fprintf(outputFileID,' Done %gs\n',toc(s2)); end;
    if verbos>1; fprintf(outputFileID,'\tStem %d of %d Main ...  Done in %gs\n',i,nStemsToTry,toc(start));   end;    
end
if verbos; fprintf(outputFileID,'Testing stems Done in %gs\n',toc(ss)); end;
% now locs contains all the motif candidates we found. 

% Calculate stats
stats=struct('mean',{},'var',{},'max',{},'min',{},'median',{},'mode',{});
nMotifs=length(locs);
for motif=1:nMotifs        
    [stats,means]=updateStatsAndMean(locs{motif},x,dPointDiffFun,stats,means,motif);
end

len=zeros(nMotifs,1); 
for ll=1:nMotifs; 
    len(ll)=length(means{ll}); 
end;
[dummy,indx]=sort(len,'descend');
locs=locs(indx); 
means=means(indx); 
stats=stats(indx);

if verbos
    tBeforeFinalization=toc(begninning);
    fprintf(outputFileID,'Total Time Before Finalization / point is: %gs\n------------------------------------------\n',tBeforeFinalization/T);    
end
end

function thismaxlen=defaultMaxLen(this,maxLen,T)
%     if this<maxLen
%         thismaxlen=maxLen;
%     else
%         thismaxlen=T;
%     end
    thismaxlen=T;
    return;
end



function [nStemsToTry,maxsmall,minlarge]=optimalDistKMeans(distLinear,optional)
    try
        idx=kmeans(distLinear(:,3),5);
    catch e
        n=size(distLinear,1);
        n1=max(1,floor(n/3)); n2=min(n-1,2*n1);
        idx=zeros(n1,1);
        idx(1:n1)=1; idx(n1+1:n2)=2;
        idx(n2+1:n)=3;
    end
    idx=idx-idx(1);
    nStemsToTry=find(idx,1,'first')-1;
    if nargout>1
        maxsmall=distLinear(nStemsToTry,3);
        idx(1:nStemsToTry)=[];
        idx=idx-idx(1);
        if nargout>2
            minlarge=distLinear(find(idx,1,'first')+nStemsToTry,3);
        end
    end
end
function [stop,isoverlapping,stemUsed,begends,meanMotif]=tryCombining(current,...
    stemLocs,stemUsed,begends,meanMotif,toreuse,minoLen,x,T,stemLen,dFun,dFunParams ...
    ,stopFun,stopFunParams,donotruese,ignoreOverlapsDuringCombination)

    stop=1;
    isoverlapping=0;
    nStemsInGroup=size(begends,1);
    for kk=1:nStemsInGroup
        if( testOverlapping([stemLocs(current),stemLocs(current)+stemLen-1] ...
                ,begends(kk,:)))
            isoverlapping=1;
            break;
        end
    end
    if(isoverlapping)
        return;
    end
    % get the next nearest stem to the one that is
    % included in the current motif
    newo=[stemLocs(current),stemLocs(current)+stemLen-1];
    % slid the stem member until it gives minimum distance                     
    % to a single stem of the mean motif    
    mind=inf; minj=0;
    for j=1:length(meanMotif)-stemLen+1                        
        dtmp=dFun(meanMotif(j:stemLen+j-1),...
                        x(newo(1):newo(2)),dFunParams);
        if(dtmp<mind)
          mind=dtmp; minj=j;
        end                        
    end                    
    newo(1)=max(1,min(T-minoLen+1,newo(1)-minj+1));
    % extend this near stem to be as large as the motif
    % we are using now
    %newo(1)=max(1,min(T-minoLen+1,newo(1)-(stemLocs(s)-begends(k,1))));
    newo(2)=newo(1)+minoLen-1;

    % slid the new location until it gives minimum 
    % diestance to the mean of the motif
    mind=inf; minj=0;
    %@todo consider a better sliding strategy like
    %continuing to slide if the best was at the edges

    %for j=[-minoLen:stemLen:-1,1:stemLen:minoLen]
    for j=[-stemLen:stemLen]
        if(newo(1)+j>0 && newo(2)+j<T+1)
            dtmp=dFun(x(newo(1)+j:newo(2)+j),meanMotif,dFunParams);
            if(dtmp<mind)
              mind=dtmp; minj=j;
            end
        else
            break;
        end
    end
    newo=newo+minj;
    % check whether the new occurrence candidate is overlapped or
    % included in an older one
    included=0; overlap=0;    
    for tmp=1:nStemsInGroup
        if newo(1)>=begends(tmp,1) && newo(2)<=begends(tmp,2)
            included=1;
            break;
        elseif ~ignoreOverlapsDuringCombination && testOverlapping(newo,begends(tmp,:))
            overlap= 1;                            
            break;
        end
    end    
    newbegends=[begends;newo];
    if ~included && ~overlap
        if(donotruese)
            [stop,toreuse]=stopFun(x,begends,newbegends,[] ...
                ,stopFunParams);                                                
        else
            [stop,toreuse]=stopFun(x,begends,newbegends,toreuse ...
                ,stopFunParams);
        end
        if stop                            
            % @here we should try to use smaller
            % regions within newo and compare them to
            % the corresponding regions in the motif.
            % if it can be fit then we can divide the
            % motif into three one includes this part of
            % the newo and the rest includes the parts
            % before it in one side and the ones after
            % it in the other.
            return;
        else            
            begends=newbegends;

            mLen=begends(:,2)-begends(:,1)+1;
            minoLen=min(mLen);  %maxoLen=max(mLen);
            meanMotif=zeros(1,minoLen);
            for kk=1:size(begends,1)
                meanMotif=meanMotif+(x(begends(kk,1):begends(kk,2)));
            end
            meanMotif=meanMotif./size(begends,1);

            stemUsed(current)=1;
            % @todo consider checking for overlaps here as
            % well            
        end
    elseif overlap                        
    end                
end
function showbegends(x,T,begends,showX,maxLen)
    for stem=1:size(begends,1)
        subplot(ceil(size(begends,1)/2),2,stem);
        hold on;            
        if showX
            range=max(1,begends(stem,1)-maxLen):min(T,begends(stem,2)+maxLen);
            plot(range,x(range),'b');
        end
        range=max(1,begends(stem,1)):min(T,begends(stem,2));
        plot(range,x(range),'r');
        drawnow;
    end
end

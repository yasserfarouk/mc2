function [locs,means,stats,maxSmallDist,minLargeDist] = ...
    sdCMD(x,lrange,candLocs,varargin)
%Shift Density Constrained Motif Discovery Algorithm
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
%           single number is passed it is treated as the maximim 
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
%                           points in the occurrences
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
% Yasser Mohammad and Toyoaki Nishida,Shift Density basec Constrained Motif
% Discovery, AAAI 2012
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%
if length(lrange)<2
    lrange=[max(3,floor(0.9*lrange)),lrange];
end

T=length(x);

if nargin<5
    chgScores=[];
end
if nargin<4
    nWindows=10; % based in results in orginial catalano's paper
end
if nargin<3
    candLocs=[];    
end


nEffectiveDists=[];   
wbar=[];
noiseLevel=.05;   % 
extraLength=0.2;
dFun=@dE;
dFunParams=[];
dPointDiffFun=@abs;         % will be used only if we need to return means and stats
windowsCreationFun=@createSWSetAround;
visualize=0;
nSkip=[];
nShifts2ReturnPerWin=1;
tryEarlyCombination=1;
minSeparation=[];
acceptableBadMatches=0.25;
nRandWinsToCalcLargeDists=max(20,min(1000,T/100));

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
            case {'nskip'}          % the number of points to skip between wbar subwindows 
                nSkip=(varargin{i+1});
            case {'nwindows','nw'} % [OPTIONAL] The number of candidate windows to use. Default is 10                                    
                nWindows=(varargin{i+1});
            case {'chgscores','chg'}     % [OPTIONAL] a time series of the same length as x. the score
                                         %             at each point is proportional to finding a motif occurrence 
                                         %             around this point. The easiest way to get this score is
                                         %             through a change point discovery routine like rsst() or cpd()
                chgScores=(varargin{i+1});            
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


%initialize other variables
w=round((1+extraLength)*lrange(2));
if isempty(wbar)
    wbar=max(3,ceil(0.1*w));
end
if isempty(nSkip)
    nSkip=max(1,floor(.5*wbar));
end
maxLen=lrange(2);
minLen=lrange(1);
minSubs=2;
nSub=numel([0:nSkip:w-wbar]);
if isempty(nEffectiveDists)
    nEffectiveDists=max(minSubs,ceil(0.5*(maxLen-wbar)./nSkip));
else
    nEffectiveDists=min(nEffectiveDists,max(minSubs,ceil((maxLen-wbar)./nSkip)));
end
if isempty(minSeparation)
    minSeparation=0;
end

% generate the noise window
%noise=x(randi(T,w,1));

% start processing
[cmpWinLocs,nWindows]=windowsCreationFun(candLocs,nWindows,w,wbar,T);


locs=[];
means=[];
stats=[];

%find nearest subsequences to noise sequences    
%[distMeansn]=compareFindBests...
%    (noise,0,x,cmpWinLocs,nWindows,nSkip,wbar,nSub,0,nEffectiveDists,dFun,dFunParams,1,1,0);
meanDistsRand=0;
for i=1:nRandWinsToCalcLargeDists
    lcs=randi(T-wbar+1,2,1);
    ddd=dFun(x(lcs(1):lcs(1)+wbar-1),x(lcs(2):lcs(2)+wbar-1),dFunParams);
    meanDistsRand=meanDistsRand+ddd;
end
meanDistsRand=meanDistsRand./nRandWinsToCalcLargeDists;
minLargeDist=meanDistsRand;         % this and the following distance are calculated per wbar points
maxSmallDist=meanDistsRand*noiseLevel;
% find nearest subsequences to candidate subsequences
ignore=zeros(nWindows,1);
for currentWin=1:nWindows       
    if ignore(currentWin)
        continue;
    end    
    %candWin,cmpWinLocs(currentWin),...
    [distMeans,bestWins,bestShifts,matchingSubs]=compareFindBests...        
        (x,cmpWinLocs,nWindows,nSkip,w,wbar,nSub,currentWin,nEffectiveDists,...
        dFun,dFunParams,minLen,nShifts2ReturnPerWin,visualize);   
    distMeans(distMeans(:,1)>maxSmallDist,:)=[];
    nAll=size(distMeans,1);
    if nAll<1
        continue;
    end    
    % remove matches with very small subs
    toremove=[zeros(1,nAll),ones(1,length(bestWins)-nAll)];
    for win=1:nAll
        if numel(matchingSubs{win,1})<minSubs
            toremove(win)=1;
        end
    end
    bestWins=bestWins(toremove~=1);
    nAll=length(bestWins);
    if nAll<1
        continue;
    end
    bestShifts=bestShifts(toremove~=1,:);
    matchingSubs=matchingSubs(toremove~=1,:);    
    distMeans=distMeans(toremove~=1,:);
    
    % matchMat will have nAll+2 rows:
    % first nAll rows will show the matching subsequence for each source
    %   subsequence    
    % row nAll+1 will have the source sub. numbers in order
    % row nAll+2 will show how many times each source sub. is matched
    matchMat=zeros(nAll+2,nSub);    
    matchMat(nAll+1,:)=1:nSub;    
    for window=1:nAll   
        list=matchingSubs{window,1}(:,1); 
        mlist=matchingSubs{window,1}(:,2); 
        nL=length(list);
        for l=1:nL
            matchMat(nAll+2,list(l))=matchMat(nAll+2,list(l))+1;
            matchMat(window,list(l))=mlist(l);
        end
    end
    matchMatSum=matchMat(nAll+2,:);
    matchMat(:,matchMatSum==0)=[];  
    matchMat(1:nAll+1,:)=(matchMat(1:nAll+1,:)-1).*nSkip;
    matchMat(matchMat==-nSkip)=NaN;
    nMatched=size(matchMat,2);        
    
    groups=zeros(nAll+1,nMatched); 
    groups(:,1)=NaN;
    groups(~isnan(matchMat(1:nAll+1,1)),1)=1;
    for win=1:nAll+1
        currentG=1;
        for m=2:nMatched
            if isnan(matchMat(win,m))
                groups(win,m)=NaN;
            else
                if matchMat(win,m)>matchMat(win,m-1)+wbar-1+minSeparation
                    currentG=currentG+1;                
                end
                groups(win,m)=currentG;
            end
        end
    end
    nG=groups(end,end);
%    matchGroups=cell(nG,1);
%    bestWindowsG=cell(nG,1);
    %bestShiftsG=cell(nG,1);    
    for g=1:nG        
        tmp=matchMat(:,groups(end,:)==g);
        gtmp=groups(:,groups(end,:)==g);
        nMyMatches=size(tmp,2);
        if nMyMatches<minSubs            
            continue;            
        end                
        % remove windows in which there is no enough matching subsequences
        isMaching=~isnan(tmp);
        badWin=sum(isMaching,2)<minSubs;
        badWin=badWin(1:end-2); % remove myseld and the counts from the end
        % remove matching windows in which the match is not contiguous        
        for win=1:nAll
            if badWin(win)
                continue;
            end
            nBad=0;                        
            nGInW=max(gtmp(win,:));
            nCont=zeros(nGInW,1);
            for m=1:nMyMatches
                if ~isMaching(win,m)
                    nBad=nBad+1;                    
                else
                    if tmp(win,m)>tmp(win,nMyMatches) || ...
                        tmp(win,m)<tmp(win,1)
                        % this match is not continuous
                        nBad=nBad+1;
                    end
                    nCont(gtmp(win,m))=nCont(gtmp(win,m))+1;
                end
            end            
            if nBad>acceptableBadMatches*nMyMatches || ...
                    max(nCont)<minSubs
                badWin(win)=1;
            end
        end
        tmp(badWin~=0,:)=[];                
        isMaching(badWin~=0,:)=[];                
        tmp=tmp(1:end-1,:);                
        nM=size(tmp,1)-1;
        if nM<1
            continue;
        end
        toDel=zeros(nMyMatches,1);
        for col=1:nMyMatches
            if any(~isMaching(:,col))
                if col==1 || toDel(col-1)
                    toDel(col)=1;                
                end
            end
        end
        for col=nMyMatches:-1:1
            if any(~isMaching(:,col))
                if col==nMyMatches || toDel(col+1)
                    toDel(col)=1;                
                end
            end
        end
        
        tmp=tmp(:,toDel~=1);
        if numel(tmp)<1
            continue;
        end
        begends=cmpWinLocs(currentWin)+[tmp(end,1),...
            tmp(end,end)+wbar-1];
        matchLen=begends(2)-begends(1);        
        bestWindowsG=bestWins(badWin==0);        
        for m=1:nM
            begends=[begends;...
                cmpWinLocs(bestWindowsG(m))+tmp(m,1),...
                cmpWinLocs(bestWindowsG(m))+tmp(m,1)+matchLen];
        end
        locs=[locs;{begends}];
        
        if ~tryEarlyCombination 
            ignore(bestWindowsG)=1;
        else
            %error('not implemented yet');
        end
        if visualize                
            figure;
            subplot(ceil((nM+1)/2),2,1);
            hold on;
            plot(x(cmpWinLocs(currentWin):cmpWinLocs(currentWin)+w-1),'g:');            
            begend=begends(1,:);
            plot(-cmpWinLocs(currentWin)+[begend(1):begend(2)],x(begend(1):begend(2))...
                ,'b','LineWidth',2);  
            for ii=1:nM
                win=bestWindowsG(ii);                
                subplot(ceil((nM+1)/2),2,ii+1);
                hold on;                            
                plot(x(cmpWinLocs(win):cmpWinLocs(win)+w-1),'r:');                
                begend=begends(ii+1,:);
                plot(-cmpWinLocs(win)+[begend(1):begend(2)],x(begend(1):begend(2))...
                    ,'r','LineWidth',2);      
            end
        end        
    end                        
end
if nargout>1
    nMotifs=length(locs);
    for i=1:nMotifs
        [stats,means]=updateStatsAndMean(locs{i},x,dPointDiffFun,stats,means,i);
    end
end
end

function [distMeans,bestWins,bestShifts,matchingSubs]=...
    compareFindBests...    
    (x,cmpLocs,nWindows,nSkip,w,wbar,nSub,currentI,nEffectiveDists,dFun,dFunParams...
    ,minLen,nShifts2ReturnPerWin,visualize)        
    
    nShifts=2*nSub-1;   %The number of possible shifts between windows
    dists=cell(nWindows-currentI,nShifts);%The distances corresponding to different shifts
    matchingSubs=cell(nWindows-currentI,nShifts);%the number of source subsequence corresponding to every distance in distances    
    %counts=zeros(nWindows,nShifts);% counts: how many times every shift is considered in dists    
    candWinStart=cmpLocs(currentI);
    for j=1:nWindows-currentI% we only consider the windows after currentI                  
        % ignore the comparison if the two windows are overlapping
        if abs(cmpLocs(currentI)-cmpLocs(currentI+j))<w
            continue;
        end
        if visualize>0 
            figure;
            allDists=zeros(nSub,nSub);
            allDists(:,:)=nan;
        end
        kk=0;
        for k=1:nSub  % the subsequences in candidate window                            
            ll=0;
            for l=1:nSub % the subsequences in comparison windows               
                cmpBeg=cmpLocs(j+currentI)+ll;
                candBeg=cmpLocs(currentI)+kk;
                if candWinStart>0 && abs(cmpBeg-candBeg)<minLen
                    ll=ll+nSkip;
                    continue;
                end
                %counts(j,l-k+nSub)=counts(j,l-k+nSub)+1;
                ddd=dFun(...
                   x(cmpBeg:cmpBeg+wbar-1),...
                   x(candBeg:candBeg+wbar-1),dFunParams);
                if visualize>0
                    allDists(k,l)=ddd;
                end
                dists{j,l-k+nSub}=[dists{j,l-k+nSub},ddd];
                matchingSubs{j,l-k+nSub}=[matchingSubs{j,l-k+nSub};k,l];                
                ll=ll+nSkip;
            end            
            kk=kk+nSkip;
        end                        
       if visualize>0           
            subplot(2,2,1);
            plot(candWinStart:candWinStart+w-1,...
                x(candWinStart:candWinStart+w-1),'LineWidth',2);
            subplot(2,2,2);
            plot(cmpLocs(j+currentI):...
                cmpLocs(j+currentI)+w-1,x(cmpLocs(j+currentI):...
                cmpLocs(j+currentI)+w-1),'LineWidth',2);
            subplot(2,2,3);
            %[X,Y]=meshgrid(1:1:nSub);
            %surf(X,Y,allDists,'FaceColor','interp',...
            %    'EdgeColor','none',...
            %    'FaceLighting','phong');                      
            sf=250/max(max(allDists));
            allDists(isnan(allDists))=255;
            image(sf.*(allDists));
            
            subplot(2,2,4);
            tmpScatter=[];
            for ss=1:nShifts
                dtmp=cell2mat(dists(j,ss))';
                ldt=length(dtmp);
                if ldt>0
                    tmpScatter=[tmpScatter;(ss-nSub).*nSkip*ones(ldt,1),dtmp];
                end
            end
            scatter(tmpScatter(:,1),tmpScatter(:,2),1);
%             dtmp=zeros(1,nShifts);
%             dtmp(:)=inf;
%             for ss=1:nShifts
%                 dtmp(ss)=mean(cell2mat(dists(j,ss)));
%             end
%             plot([-nSub+1:nSub-1],dtmp,'LineWidth',2);
        end
    end      
    
    % now we have a list of distances at every window for every shift. 
    % we would like to find the best windows and the best shift for each.
    % we have also to keep the beginning of the matching segment in the
    % source window through begends
    distMeans=zeros(nWindows-currentI,nShifts);  
    bestShifts=zeros(nWindows-currentI,nShifts);    
    for j=1:nWindows-currentI
        for k=1:nShifts            
            if length(dists{j,k})<1
                distMeans(j,k)=inf;
                continue;
            end
            [tmp,I]=sort(squeeze(dists{j,k}),'ascend');            
            matchingSubs{j,k}=matchingSubs{j,k}(I,:);            
            max2Consider=min(length(tmp),nEffectiveDists);
            matchingSubs{j,k}=matchingSubs{j,k}(1:max2Consider,:);            
            distMeans(j,k)=sum(tmp(1:max2Consider))./max2Consider;
            %dists{j,k}=dists{j,k}(I(1:max2Consider));
        end
        [distMeans(j,:),bestShifts(j,:)]=sort(distMeans(j,:),'ascend');
        matchingSubs(j,:)=matchingSubs(j,bestShifts(j,:));        
    end    
    [distMeans,bestWins]=sortrows(distMeans,[1:nShifts]);
    if nargout>1
        bestShifts=(bestShifts(bestWins,:)-nSub).*nSkip;
        matchingSubs=matchingSubs(bestWins,:);        
        bestWins=bestWins+currentI;    
        if nargin>12
            bestShifts=bestShifts(:,nShifts2ReturnPerWin);
            distMeans=distMeans(:,nShifts2ReturnPerWin);
            matchingSubs=matchingSubs(:,nShifts2ReturnPerWin);            
            nShiftsMax=nShifts2ReturnPerWin;
        else
            nShiftsMax=nShifts;
        end
        for j=1:nWindows-currentI
            for k=1:nShiftsMax
                if ~isempty(matchingSubs{j,k})
                    [matchingSubs{j,k},I]=sortrows(matchingSubs{j,k},1);                
                end
            end
        end
    end
end


function [loc,nWindows]=createSWSetCenter(candLoc,nWindows,w,wbar,T)
% creates and returns a subwindow set in the form of a nWindows*nSub*2 matrix    
    loc=zeros(nWindows,1);
    for i=1:nWindows
        loc(i)=min(T-w+1,max(1,candLoc(i)-floor(w/2)));                    
    end
end
function [loc,nWindows]=createSWSetAround(candLoc,nWindows,w,wbar,T)
% creates and returns a subwindow set in the form of a (nWindows)*nSub*2 matrix
% notice that nWindows in the output can be upt to 2*nWindows [input]
    loc=zeros(2*nWindows,1);
    for i=1:nWindows
        loc(2*(i-1)+1)=max(1,min(T-w+1,candLoc(i)-w));
        loc(2*i)=max(1,min(T-w+1,candLoc(i)+wbar));
    end
    loc=unique(loc);
    nWindows=length(loc);
end
% function d=dE(s1,s2,optional)
% % Distance with possible normlaization. it automatically compresses/extends the
% % second series to be the length of the first. 
%     normalize=1;
%     power=2;
%     mn=min(length(s2),length(s1));
%     s1=compressSeries(s1,mn);%,max(s1),min(s1));
%     s2=compressSeries(s2,mn);%,max(s1),min(s1));
%     if nargin>2 && ~isempty(optional)
%         normalize=optional{1};
%         if length(optional)>1
%             power=optional{2};
%         end
%     end    
%     if normalize
%         %if(sum(s1.^2)>eps)        
%             %s1=s1./sqrt(sum(s1.^2));
%             delta=max(s1)-min(s1);
%             if delta>0.000000000001
%                 s1=s1./(delta);
%             end
%         %end
%         %if(sum(s2.^2)>eps)
%             delta=max(s2)-min(s2);
%             if delta>0.000000000001
%                 s2=s2./(delta);
%             end
%             %s2=s2./sqrt(sum(s2.^2));
%         %end
%     end
%     d=mean(abs((s2-s1).^power));
% end

% function destMachingSubs=findMatchingSubs(bestShifts,matchingSubs,nAll)
%     destMachingSubs=matchingSubs;
%     nShifts=size(bestShifts,2);
%     for shift=1:nShifts
%         for i=1:nAll
%             destMatchingSubs{i,shift}=destMatchingSubs{i,shift}+...
%                 bestShifts(i,shift);
%         end
%     end
% end
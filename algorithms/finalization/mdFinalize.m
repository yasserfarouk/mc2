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

function [locs,means,stats] = mdFinalize(x,lrange,candLocs,locs,means,stats,varargin)
%Finalizes the motif discovery operation by combining, extending and
%further processing of motifs
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
%           the changes. If empty then it will be ignored
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
% varargin  See the switch statement for their meaning
%
% Outputs:
% ========
% same as inputs with the same names but after finalization steps
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

% process inputs
if length(lrange)<2
    lrange=[max(3,floor(0.9*lrange)),lrange];
end

% intialize internal parameters
stemLen=[];
maxLen=[];
minLen=[];
epsilon=[];
acceptableMeanIncrease=[];

visualize=0;
outputFileID=1;
verbos=0;
maxFinalizationTrials=5;
maxSlidingRetrialsPerExtension=3;
minExtension=[];             % minimum number of points to extend 

useBisection=1;
tryRemoveExtraStems=1;
trySlidingDuringFinalExtension=1;
tryFullDetectionDuringFinalization=1;
tryCombiningSimilarMeanMotifsDuringFinalization=1;
tryRemovingDuplicatesDuringFinalization=1;

maxSmallDist=inf;
minLargeDist=inf;
meanDistFraction4Combine=0.5;
maxMeanErrsToTrySliding=1.1;
maxOverlapWithinMotif=0;
maxOverlapBetweenMotifs=1-.00000001;
significance=0.05;
removeAllZeros=1;
removeAllConstant=1;
ignoreOverlapsDuringExtension=1;
ignoreOverlapsDuringFinalization=1;
ignoreOverlappingShorterMotifsDuringDetection=1;
useCandLocsDuringDetection=~isempty(candLocs);
minNOccurrences=2;
maxNOccurences=10;

dPointDiffFun=@abs;
dFun=@dE;
dFunParams=[];
stopFun=@stopIfAboveMean;
stopFunParams=[];
stopFunParamsForFinalExtension=[];
stopFunParamsStrict=[];             % only accepts a change if the new points will have virtually no negative effect at all
maxLenFun=@defaultMaxLen;
motifDetectionFun=@detectMotifDefault;
motifDetectionFunParams=[];
% change internal parameters using optional arguments

nArgs=size(varargin,2);
if(nArgs==1 && iscell(varargin{1}))
    varargin=varargin{1};
    nArgs=size(varargin,2);
end
if(nArgs>0)
    if(mod(nArgs,2)~=0)
        error('The optional arguments must be in the form name,value so they must be even!!!');
    end
    for i=1:2:nArgs
        switch(lower(varargin{i}))                              
            case {'minlargedist'}          
                minLargeDist=(varargin{i+1});                   
            case {'maxsmalldist'}          
                maxSmallDist=(varargin{i+1});                   
            case {'tryremoveextrastems'}          
                tryRemoveExtraStems=(varargin{i+1});                   
            case {'tryslidingduringfinalextension'}          
                trySlidingDuringFinalExtension=(varargin{i+1});                   
            case {'tryfulldetectionduringfinalization'}          
                tryFullDetectionDuringFinalization=(varargin{i+1});                   
            case {'trycombiningsimilarMeanmotifsduringfinalization'}          
                tryCombiningSimilarMeanMotifsDuringFinalization=(varargin{i+1});                   
            case {'tryremovingduplicatesduringfinalization'}          
                tryRemovingDuplicatesDuringFinalization=(varargin{i+1});                               
            case {'usebisection','bisect'}          % whether or not to use bisection when trying to fastly extend stems
                useBisection=(varargin{i+1});                   
            case {'ignoreoverlapsduringextension'}          
                ignoreOverlapsDuringExtension=(varargin{i+1});                               
            case {'ignoreoverlapsduringfinalization'}          
                ignoreOverlapsDuringFinalization=(varargin{i+1});                   
            case {'ignoreoverlappingshorterMotifsduringdetection'}          
                ignoreOverlappingShorterMotifsDuringDetection=(varargin{i+1});                   
            case {'usecandlocsduringdetection'}          
                useCandLocsDuringDetection=(varargin{i+1});                       
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
            case {'stopfunparamsforfinalextension'}          
                stopFunParamsForFinalExtension=(varargin{i+1});                   
            case {'stopfunparamsstrict'}          
                stopFunParamsStrict=(varargin{i+1});                   
            case {'maxlenfun'}          
                maxLenFun=(varargin{i+1});                   
            case {'motifdetectionfun'}          
                motifDetectionFun=(varargin{i+1});                   
            case {'motifdetectionfunparams'}          
                motifDetectionFunParams=(varargin{i+1});                                        
            case {'removeallzeros'}          
                removeAllZeros=(varargin{i+1});                   
            case {'removeallconstant'}          
                removeAllConstant=(varargin{i+1});                               
            case {'minnoccurrences'}          
                minNOccurrences=(varargin{i+1});                   
            case {'maxnoccurences'}          
                maxNOccurences=(varargin{i+1});            
            case {'acceptablemeanincrease'}          
                acceptableMeanIncrease=(varargin{i+1});                   
            case {'visualize'}          
                visualize=(varargin{i+1});                   
            case {'outputfileid'}          
                outputFileID=(varargin{i+1});                   
            case {'verbos'}          
                verbos=(varargin{i+1});                   
            case {'maxfinalizationtrials'}          
                maxFinalizationTrials=(varargin{i+1});                   
            case {'maxslidingretrialsperextension'}          
                maxSlidingRetrialsPerExtension=(varargin{i+1});                   
            case {'minextension'}          
                minExtension=(varargin{i+1});              
            case {'useminmaxdists'}          
                useMinMaxDists=(varargin{i+1});                   
            case {'meandistfraction4combine'}          
                meanDistFraction4Combine=(varargin{i+1});                   
            case {'maxMeanErrsToTrySliding'}          
                maxMeanErrsToTrySliding=(varargin{i+1});                   
            case {'maxoverlapwithinmotif'}          
                maxOverlapWithinMotif=(varargin{i+1});                   
            case {'maxoverlapbetweenmotifs'}          
                maxOverlapBetweenMotifs=(varargin{i+1});                   
            case {'significance'}          
                significance=(varargin{i+1});                               
            case {'stemlen','sl'}          
                stemLen=(varargin{i+1});                                                   
            case {'epsilon'}          
                epsilon=(varargin{i+1});                               
            case {'maxlen'}          
                maxLen=(varargin{i+1});                       
            case {'minlen'}          
                minLen=(varargin{i+1});                               
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
if isempty(stemLen)
    stemLen=floor(lrange(1)/2);
end
if isempty(maxLen)
    maxLen=lrange(2);
end
if isempty(minLen)
    minLen=lrange(1);
end
if isempty(stopFunParams)
    stopFunParams={significance,epsilon,dPointDiffFun,acceptableMeanIncrease};
end
if isempty(stopFunParamsForFinalExtension)
    stopFunParamsForFinalExtension={significance,epsilon,dPointDiffFun,0.5*acceptableMeanIncrease};
end
if isempty(stopFunParamsStrict)
    stopFunParamsStrict={significance/10,epsilon,dPointDiffFun,acceptableMeanIncrease/100};
end
if isempty(minExtension)
    minExtension=1;
end
useMinMaxDists=~isinf(maxSmallDist);
if(iscell(stats))
    stats=cell2mat(stats);
end
if isempty(means)
    nMotifs=length(locs);
    for i=1:nMotifs
        [stats,means]=updateStatsAndMean(locs{i},x,dPointDiffFun,stats,means,i);
    end
end

nFinalizationTrials=0;
somethingMayHaveChanged=100*ones(length(locs),1);
removedMotifs=[];
% 1 = combined with other motif
% 2 = some occurrences were removed for overlapping (may be some where
%       added)
% 4 = Some occurrences were slided
% 8 = The motif was extended
% 16 = Some new occurrences found during detection
% 100 = initialization value
while nFinalizationTrials<maxFinalizationTrials && ...
        (any(somethingMayHaveChanged~=0) || ~isempty(removedMotifs))
    if visualize
        figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        colors={'r','b','k','c','m','y','--r','--b','--g','--k','--c','--m'...
            ,'-.r','-.b','-.y','-.k','-.c','-.m',':r',':b',':y',':k',':c',':m'};
        nMotifs=length(means);        
        nxt=1;
        nRows=ceil(nMotifs/2)+1;
        mx=max(max(x))+.05.*max(max((x))); 
        mn=min(min(x))-.05.*min(min(x));
        for m=1:nMotifs
            subplot(nRows,2,nxt);
            hold on;
            no=size(locs{m},1);
            nMax=max([numel(means{m}),(max(locs{m}(:,2)-locs{m}(:,1)+1))]);
            for i=1:no
                plot(compressSeries(x(locs{m}(i,1):locs{m}(i,2)),nMax),':g');
            end
            plot(compressSeries(means{m},nMax),colors{mod(m,length(colors))});
            ylim([mn mx]);
            title(sprintf('N(%d),Length (%d,%d),mean(%g) var(%g) \nmax(%g) median(%g) mode(%g) min(%g)' ...
                ,size(locs{m},1),(min(locs{m}(:,2)-locs{m}(:,1)+1)),(max(locs{m}(:,2)-locs{m}(:,1)+1))...
                ,stats(m).mean,stats(m).var,stats(m).max...
                ,stats(m).median,stats(m).mode,stats(m).min));
            nxt=nxt+1;
        end
        if mod(nxt,2)<0.5
            nxt=nxt+1;
        end
        subplot(nRows,2,[nxt:nxt+1]);
        locv=zeros(size(x));
        locv(candLocs)=max(max(x));
        locv(locv<0.5)=inf;
        plot(x,':g');  hold on; plot(locv,'*r'); hold off;
        nMotifs=length(means);

        hold on;
        for i=1:nMotifs    
            tmp=zeros(size(x)); tmp(tmp<1)=inf;
            clr=colors{mod(i,length(colors))};
            for j=1:size(locs{i},1);
                tmp(locs{i}(j,1):locs{i}(j,2))=x(locs{i}(j,1):locs{i}(j,2));
            end
            plot(tmp,clr);
        end
        xlabel(sprintf('Finalization trial %d of %d', ...
            nFinalizationTrials,maxFinalizationTrials));
    end
    if verbos; 
        fprintf(outputFileID,'Finalization steps [%d of %d] Trials ...:\n',...
            nFinalizationTrials+1,maxFinalizationTrials);
        if nFinalizationTrials>0
            fprintf(outputFileID,'\tChanged (of %d):',length(locs));
            for i=1:length(locs)
                if somethingMayHaveChanged(i)
                    fprintf(outputFileID,'%d[%d,%d] ', i,...
                        somethingMayHaveChanged(i),size(locs{i},1));
                end
            end
            fprintf(outputFileID,'\n');
        end
        start1=tic; 
    end;    
    removedMotifs=[];
    recentChanges=somethingMayHaveChanged;
    somethingMayHaveChanged=zeros(length(locs),1);    
    nFinalizationTrials=nFinalizationTrials+1;
    if verbos; 
        fprintf(outputFileID,'\tCombining Motifs, removing duplicates and sigelton motifs [%d of %d Trials]... ' ...
            ,nFinalizationTrials,maxFinalizationTrials); 
        start=tic; 
    end;    
    nMotifs=length(locs);
    % Combine all motifs that have very similar means:
    %     Find motifs that need to be combined (having means that are
    %     similar both using pointdiff function and distance function)
    if tryCombiningSimilarMeanMotifsDuringFinalization
        tocombine=[];
        for motif=1:nMotifs-1        
            meandist=stats(motif).mean; vardist=stats(motif).var;
            mean1=means{motif}; nm1=length(mean1);
            for motif2=motif+1:nMotifs
                if ~recentChanges(motif) && ~recentChanges(motif2)
                    continue;
                end
                if length(means{motif2})~=nm1
                    continue;
                end
                mean2=(means{motif2});
                diffs=dPointDiffFun(mean2-mean1);
                h1=safeztest(diffs,meandist,vardist,significance,'right',nm1*epsilon);
                h2=safeztest(diffs,stats(motif2).mean,stats(motif2).var...
                    ,significance,'right',nm1*epsilon);
                if useMinMaxDists
                    meanDSmall=dFun(means{motif},means{motif2},dFunParams)...
                        <meanDistFraction4Combine*maxSmallDist;
                else
                    meanDSmall=1;
                end
                if ~h1 && ~h2 && meanDSmall
                    nMotifs2=size(locs{motif2},1);                    
                    l1=size(locs{motif},1); l2=size(locs{motif2},1);
                    mindis=inf;
                    for k=1:l1
                        for kk=1:l2
                            ddd=abs(locs{motif}(k,1)-locs{motif2}(kk,1));
                            if ddd<mindis; mindis=ddd; end;
                        end
                    end
                    if mindis>minLen
                        tocombine=[tocombine;motif,motif2];                
                    end
                end
            end
        end
        %     Combine Them
        toremovem=[];
        torecalcstats=[];
        if ~isempty(tocombine)        
            for i=1:size(tocombine,1)
                somethingMayHaveChanged(tocombine(i,1))=1;
                somethingMayHaveChanged(tocombine(i,2))=1;
                recentChanges(tocombine(i,1))=1;
                recentChanges(tocombine(i,2))=1;
                n1=length({tocombine(i,1)}); n2=length({tocombine(i,2)});
                means{tocombine(i,1)}=(means{tocombine(i,1)}.*n1+ ...
                    (means{tocombine(i,2)}).*n2)./(n1+n2);
                locs{tocombine(i,1)}=[locs{tocombine(i,1)};locs{tocombine(i,2)}];    
                toremovem=[toremovem;tocombine(i,2)];            
                torecalcstats=[torecalcstats;tocombine(i,1)];
            end
            for i=motif:length(torecalcstats)
                [stats,means]=updateStatsAndMean(locs{motif},...
                    x,dPointDiffFun,stats,means,motif);
            end
            means(toremovem)=[];
            locs(toremovem)=[];
            stats(toremovem)=[];
            somethingMayHaveChanged(toremovem)=[];
            recentChanges(toremovem)=[];
            removedMotifs=[removedMotifs,toremovem(:)'];
        end
        % Combine all motifs that have very similar means Done
    end
    if tryRemovingDuplicatesDuringFinalization
        % Remove duplicates within each motif. To do that we remove any
        % occurrence that is exactly the same as another or that is overlapping
        % another which is nearer to the mean
        nMotifs=length(locs);
        for motif=1:nMotifs
            if ~recentChanges(motif)
                continue;
            end
            begends=locs{motif};
            n=size(begends,1);
            toremove=zeros(n,1);
            for i=1:n-1
                for j=i+1:n
                    if(~toremove(j))
                        if begends(i,1)==begends(j,1) && begends(i,2)==begends(j,2)
                            toremove(j)=1;
                        elseif  testOverlapping(begends(i,:),begends(j,:))
                            %overlap 
                            ovelapfraction=overlapfraction(begends(i,:),begends(j,:));
                            if ovelapfraction>maxOverlapWithinMotif
                                di=dFun(x(begends(i,1):begends(i,2)),means{motif},dFunParams);
                                dj=dFun(x(begends(j,1):begends(j,2)),means{motif},dFunParams);
                                if(dj<di)
                                    toremove(i)=1;
                                else
                                    toremove(j)=1;
                                end
                            end
                        end
                    end
                end
            end
             if any(toremove)
                 recentChanges(motif)=2;
             end
            begends(toremove~=0,:)=[];
            [stats,means]=updateStatsAndMean(locs{motif},...
                    x,dPointDiffFun,stats,means,motif);
            [dummy,indx]=sort(begends(:,1));    
            locs{motif}=begends(indx,:);        
        end
        % Remove duplicates within each motif DONE        
    end
    
    % remove motifs with a single occurrence
    toremove=zeros(nMotifs,1);
    for motif=nMotifs:-1:1
        if ~recentChanges(motif)
            continue;
        end
        if(size(locs{motif},1)<2)
            toremove(motif)=1;
        end
    end
    means(toremove~=0)=[];
    locs(toremove~=0)=[];
    stats(toremove~=0)=[];
    somethingMayHaveChanged(toremove~=0)=[];
    recentChanges(toremove~=0)=[];
    removedMotifs=[removedMotifs,find(toremove~=0)'];
    
    % remove motifs with a single occurrence DONE
    if verbos; fprintf(outputFileID,'Done in %gs\n',toc(start)); end;    

    if verbos>1; fprintf(outputFileID,'\tCalculating stats and sorting motifs ...'); start=tic; end;                
    stats=struct('mean',{},'var',{},'max',{},'min',{},'median',{},'mode',{});
    nMotifs=length(locs);
    for motif=1:nMotifs
        begends=locs{motif};
        [stats,means]=updateStatsAndMean(begends,x,dPointDiffFun,stats,means,motif);
    end

    % sort the motifs according to mean distance
    %meandists=[stats.mean];
    %[dummy,indx]=sort(meandists,'ascend');
    %locs=locs(indx); 
    %means=means(indx); 
    %stats=stats(indx);
    % sort the motifs according to length
    len=zeros(nMotifs,1); 
    for ll=1:nMotifs; 
        len(ll)=length(means{ll}); 
    end;
    [dummy,indx]=sort(len,'descend');
    locs=locs(indx); 
    means=means(indx); 
    stats=stats(indx);
    if verbos>1; fprintf(outputFileID,'Done in %gs\n', toc(start)); end;
    if verbos>1; fprintf(outputFileID,'\tResolving overlaps between different motifs ...'); start=tic; end;
    
    % check for overlaps betwen occurrences of different motifs and resovle
    % them. Notice that always the motif with lower index is 'better' and
    % should never be removed
    nMotifs=length(locs);
    toremove=zeros(nMotifs,1);    
    for i=1:nMotifs-1
        if ~recentChanges(i)
            continue;
        end
        ni=size(locs{i},1);
        meanMotif=means{i};
        begends=locs{i};
        for j=i+1:nMotifs
            if(toremove(j))
                continue;
            end
            nj=size(locs{j},1);     % remember that within each motif occurrences are ordered
            overlapwith=zeros(nj,2);
            for oj=1:nj
                for oi=1:ni;
                    if(locs{i}(oi,1)>locs{j}(oj,2))
                        break;
                    elseif locs{i}(oi,1)==locs{j}(oj,1) && locs{i}(oi,2)==locs{j}(oj,2)
                        overlapwith(oj,1)=oi; overlapwith(oj,2)=1;
                    elseif testOverlapping(locs{i}(oi,:),locs{j}(oj,:))
                        overlapwith(oj,1)=oi;
                        overlapwith(oj,2)=overlapfraction(locs{i}(oi,:),locs{j}(oj,:));
                    end
                end
            end
            % if we found an overlap we will remove this later motif (j) but after
            % trying to combine its nonoverlapping members with the motif i
            toremove(j)=mean(overlapwith(:,2))>=maxOverlapBetweenMotifs;
            if(toremove(j))
                for oj=1:nj
                    if(overlapwith(oj,1)<0.5)                
                        % try without changing the length
                        trialLen=floor(max([stemLen,0.5*(locs{j}(oj,2)-locs{j}(oj,1)+1)]));
                        beg=floor(0.5*(locs{j}(oj,1)+locs{j}(oj,2))-0.5*trialLen);
                        newo=[beg,beg+trialLen-1];
                        % slid the stem member until it gives minimum distance                     
                        % to a single stem of the mean motif    
                        mind=inf; minj=0;
                        for s=1:length(meanMotif)-trialLen+1                        
                            dtmp=dFun(meanMotif(s:trialLen+s-1),...
                                            x(newo(1):newo(2)),dFunParams);
                            if(dtmp<mind)
                              mind=dtmp; minj=j;
                            end                        
                        end                    
                        newo(1)=max(1,min(T-numel(meanMotif)+1,newo(1)-minj+1));
                        % extend this near stem to be as large as the motif
                        % we are using now                
                        newo(2)=newo(1)+numel(meanMotif)-1;
                        newbegends=[begends;newo];                                    
                        stop=stopFun(x,begends,newbegends,[],stopFunParams);                                                                
                        if stop
                        else
                            begends=newbegends;
                            [stats,means]=updateStatsAndMean...
                                (begends,x,dPointDiffFun,stats,means,i);
                            locs{i}=begends;
                        end
                    end
                end 
            end
        end
    end    
    means(toremove~=0)=[];
    locs(toremove~=0)=[];
    stats(toremove~=0)=[];
    somethingMayHaveChanged(toremove~=0)=[];
    recentChanges(toremove~=0)=[];
    removedMotifs=[removedMotifs,find(toremove~=0)'];
    someMotifsRemovedBecauseOverlapping=any(toremove);
    if verbos>1; fprintf(outputFileID,'Done in %gs[removed %d motifs]\n', toc(start),sum(toremove)); end;
    % removing unnecessary motifs overlapping with others ... DONE
    
    % trying to adjust the motif by removing unnecessary outer parts and
    % then extending it to the maximum possible
    if verbos>1; fprintf(outputFileID,'\tRemoving Extra stems and possibly extending motifs ...'); start=tic; end;
    if tryRemoveExtraStems && (nFinalizationTrials==1 || someMotifsRemovedBecauseOverlapping)
        nMotifs=length(locs);
        % @todo it will be faster if rather than extending from the middle
        % we try to remove one stem at the time from the ends but this will
        % require us to write another stopFun and may be modify extendStem
        % because currently stopFun assumes that either a new begend is
        % added or an old one is extended. we will need another one like
        % continueRemovalFun or so to replace stopFun
        for motif=1:nMotifs
            if ~recentChanges(motif)
                continue;
            end
            begends=locs{motif};            
            n=size(begends,1);
            if verbos>2; fprintf(outputFileID,'\t\tSliding occurrences to best position[%d of %d] ....',motif,nMotifs); start=tic; end;
            tryaSlidingAgain=trySlidingDuringFinalExtension && ...
                stats(motif).mean<=maxMeanErrsToTrySliding*stats(1).mean && ...
                (~useMinMaxDists || ...
                 ~isMeanDistBetweenOccurencesSmall(x,locs{motif},maxSmallDist,minLargeDist,dFun,dFunParams));
            nSlidingTrialsPerExtension=0;            
            while tryaSlidingAgain && ...
                    (nSlidingTrialsPerExtension<maxSlidingRetrialsPerExtension)
                tryaSlidingAgain=0;
                nSlidingTrialsPerExtension=nSlidingTrialsPerExtension+1;
                if verbos>2; fprintf(outputFileID,'Trial:%d of %d ',nSlidingTrialsPerExtension,maxSlidingRetrialsPerExtension); end;
                for i=1:n
                    IwasMovedLeft=0; myOriginalBegend=begends(i,:);
                    newbegends=begends;
                    % Try to move one motif to the left as long
                    % as there is no increase in the mean error
                    newbegends(i,:)=newbegends(i,:)-stemLen;
                    if(newbegends(i,1)<1)
                        continue;
                    end
                    [stop,toreuse]=stopFun(x,begends,newbegends,[],stopFunParamsStrict);
                    while(~stop)                                       
                        begends(i,:)=newbegends(i,:);
                        IwasMovedLeft=1;
                        % check that no occurrences overlap
                        overlap=0; 
                        for k=[1:i-1,i+1:n]                            
                            if ~ignoreOverlapsDuringFinalization && testOverlapping(begends(i,:),begends(k,:))
                                overlap=1; break;
                            end
                        end                        
                        % if there is an overlap, we just stop and exit
                        if overlap
                            % @todo do something better in case of overlap!!!!!!!
                            %warning('DGR:OVERPLAP','Overlap detected ... DGR');
                            break;
                        end
                        newbegends(i,:)=newbegends(i,:)-stemLen;
                        if(newbegends(i,1)<1)
                            break;
                        end                        
                        [stop,toreuse]=stopFun(x,[],newbegends,toreuse,stopFunParamsStrict);
                    end
                    myBestLeft=begends(i,:);                    
                    
                    % try to move to the right from original location
                    IwasMovedRight=0;
                    newbegends=begends;
                    newbegends(i,:)=myOriginalBegend+stemLen;
                    if(newbegends(i,2)>T)
                        continue;
                    end
                    [stop,toreuse]=stopFun(x,[],newbegends,toreuse,stopFunParamsStrict);
                    while(~stop)
                        begends(i,:)=newbegends(i,:);
                        IwasMovedRight=1;                        
                        % check that no occurrences overlap
                        overlap=0;                        
                        for k=[1:i-1,i+1:n]                            
                            if ~ignoreOverlapsDuringFinalization && testOverlapping(begends(i,:),begends(k,:))
                                overlap=1; break;
                            end
                        end                        
                        % if there is an overlap, we just stop and exit
                        if overlap
                            % @todo do something better in case of overlap!!!!!!!
                            %warning('DGR:OVERPLAP','Overlap detected ... DGR');
                            break;
                        end
                        newbegends(i,:)=newbegends(i,:)+stemLen;
                        if(newbegends(i,2)>T)
                            break;
                        end                        
                        [stop,toreuse]=stopFun(x,[],newbegends,toreuse,stopFunParamsStrict);
                    end
                    myBestRight=begends(i,:);                    
                    if IwasMovedLeft && IwasMovedRight
                        dLeft=dFun(x(begends(i,1):begends(i,2)),...
                                   x(myBestLeft(1):myBestLeft(2)), dFunParams);
                        dRight=dFun(x(begends(i,1):begends(i,2)),...
                                   x(myBestLeft(1):myBestLeft(2)), dFunParams);
                        if dLeft<dRight
                            begends(i,:)=myBestLeft;
                        else
                            begends(i,:)=myBestRight;
                        end
                    elseif IwasMovedLeft
                        begends(i,:)=myBestLeft;
                    elseif IwasMovedRight
                        begends(i,:)=myBestRight;
                    end
                    if IwasMovedLeft || IwasMovedRight                        
                        tryaSlidingAgain=1;
                        locs{motif}=begends;
                    end
                end 
                if tryaSlidingAgain
                    somethingMayHaveChanged(motif)=somethingMayHaveChanged(motif)+4;
                    recentChanges(motif)=4;
                end
            end
            if verbos>2; fprintf(outputFileID,'Done in %gs\n', toc(start)); end;
            if verbos>2; fprintf(outputFileID,'\t\tExtending from the middle out [%d of %d] ....',motif,nMotifs); start=tic; end;
            myStemLen=min(stemLen,min(begends(:,2)-begends(:,1)+1));
            thismaxLen=maxLenFun(max(begends(:,2)-begends(:,1)+1),maxLen,T);
            centers=floor(0.5*(begends(:,2)+begends(:,1)));
            begends(:,1)=max(round(centers-0.5*myStemLen),1);
            begends(:,2)=begends(:,1)+myStemLen-1;
            [begends,toreuse]=extendStem([-myStemLen,myStemLen],x,T,thismaxLen,begends ...
                ,stopFun,[],stopFunParamsForFinalExtension,maxSmallDist,minLargeDist...
                ,dFun,dFunParams,minExtension,1,useBisection,0,ignoreOverlapsDuringExtension);
            if verbos>2; fprintf(outputFileID,'Done in %gs\n', toc(start)); end;            
            chng=((max(max(abs(begends-locs{motif})))>0.5));
            if chng 
                somethingMayHaveChanged(motif)=somethingMayHaveChanged(motif)+8;
                recentChanges(motif)=8;
                locs{motif}=begends;                           
                [stats,means]=updateStatsAndMean...
                    (begends,x,dPointDiffFun,stats,means,motif);
            end
        end        
    end
    if verbos>1; fprintf(outputFileID,'Done in %gs\n', toc(start)); end;    
    
    % try detection of motifs in the complete set
    if verbos>1; fprintf(outputFileID,'\tFull series detection ...:\n'); start=tic; end;
    if tryFullDetectionDuringFinalization
        % order motifs with length. The longer the motif the better is it
        % for us .... may be we can also consider the statistics of the
        % motif but the rationale here is that the motif would not have
        % been allowed to become long if its statistics were bad
        len=zeros(nMotifs,1); for ll=1:nMotifs; len(ll)=length(means{ll}); end;
        [dummy,indx]=sort(len,'descend');
        locs=locs(indx); 
        means=means(indx); 
        stats=stats(indx);
        isChanged=0;
        toremove=zeros(nMotifs,1);
        for i=1:nMotifs
            if ~recentChanges(i)
                continue;
            end
            if verbos>2; fprintf(outputFileID,'\t\tFull series detection [%d of %d] ....',i,nMotifs); start=tic; end;
            % check that the motif occurrences are not included in any of
            % the previously checked motifs            
            if ~ignoreOverlappingShorterMotifsDuringDetection
                no=size(locs{i},1);
                included=zeros(no,1);
                for j=1:i-1
                    for k=1:no        
                        if included(k) 
                            continue;
                        end
                        no2=size(locs{j},1);
                        for k2=1:no2
                            over=testOverlapping(locs{i}(k,:),locs{j}(k2,:));
                            if over 
                                included(k)=1;
                                break;
                            end
                        end
                    end
                end
                % if any occurrence is overlapping with a better motif then do
                % not check this motif out
                if any(included)
                    if verbos>2; fprintf(outputFileID,'Continued in %gs\n', toc(start)); end;    
                    continue;
                end
            end
            % now we know this motif is genuine and we would like to find
            % it all over the time series
            if size(locs{i},1)>maxNOccurences
                toremove(i)=1;                
            else
                if useCandLocsDuringDetection
                    [locs{i},isChanged]=motifDetectionFun(x,locs{i},means{i},...
                        stats(i),candLocs,maxSmallDist,minLargeDist,dFun,dFunParams...
                            ,stopFun,stopFunParamsStrict,motifDetectionFunParams);            
                else
                    [locs{i},isChanged]=motifDetectionFun(x,locs{i},means{i},...
                        stats(i),[],maxSmallDist,minLargeDist,dFun,dFunParams...
                            ,stopFun,stopFunParamsStrict,motifDetectionFunParams);            
                end
                if size(locs{i},1)>maxNOccurences || size(locs{i},1)<minNOccurrences
                    toremove(i)=1;
                else
                    if isChanged
                        somethingMayHaveChanged(i)=somethingMayHaveChanged(i)+16;
                        recentChanges(i)=16;
                        [stats,means]=updateStatsAndMean(locs{i},x,dPointDiffFun...
                            ,stats,means,i);
                    end 
                end
            end
            if verbos>2; fprintf(outputFileID,'Done in %gs [%d,%d]\n',...
                    toc(start),isChanged,size(locs{i},1)); end;    
        end
        
        % remove ubiqutous and very rare motifs
        means(toremove~=0)=[];
        locs(toremove~=0)=[];
        stats(toremove~=0)=[];
        somethingMayHaveChanged(toremove~=0)=[];
        recentChanges(toremove~=0)=[];
        removedMotifs=[removedMotifs,find(toremove~=0)'];   
    end
    if verbos; fprintf(outputFileID,'Finalizing motifs and removing unnecessary ends Done in %gs [nFinalizationTrials=%d]\n',toc(start1)...
            ,nFinalizationTrials); end;
end

% Calculate stats again (should be unnecessary)
if nargout>1
    stats=struct('mean',{},'var',{},'max',{},'min',{},'median',{},'mode',{});
    nMotifs=length(locs);
    for motif=1:nMotifs        
        [stats,means]=updateStatsAndMean(locs{motif},x,dPointDiffFun,stats,means,motif);
    end
end
% remove all-zero motifs if needed
if removeAllZeros
   nMotifs=length(locs);
   toremove=[];
   for i=1:nMotifs
       if sum(means{i}.^2)<1e-9*length(means{i})
           toremove=[toremove;i];
       end
   end
   locs(toremove)=[];
   means(toremove)=[];
   stats(toremove)=[];
end
if removeAllConstant
   nMotifs=length(locs);
   toremove=[];
   for i=1:nMotifs
       if sum((means{i}(2:end)-means{i}(1:end-1)).^2)<1e-9*(length(means{i}-1))
           toremove=[toremove;i];
       end
   end
   locs(toremove)=[];
   means(toremove)=[];
   stats(toremove)=[];
end
if verbos; 
    tEnd=toc(begninning);
    fprintf(outputFileID,'Total Finalization time / point %gs\n',tEnd/T);
    fprintf(outputFileID,'Total time / point %gs\n',(tEnd+tBeforeFinalization)/T); 
end;
end

function thismaxlen=defaultMaxLen(this,maxLen,T)
    thismaxlen=T;
    return;
end


function h=safeztest(diffs,m,v,significance,type,epsilon)
    newm=median(diffs);
    if newm<m
        h=0;
    else
        if abs(v)<epsilon
            h=(newm-m>epsilon);
        else
            h=ztest(diffs,m,v,significance,type);
            if isnan(h)
                warning('DGR:ZTESTNAN','z-test failed. we will just compare means');
                h=(newm-m>epsilon);
            end            
        end
    end
end




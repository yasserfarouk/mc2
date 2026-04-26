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

% Tests differnt  motif discovery algorithms
% See testcmdLNoises.m for the parameters that need to be set before using
% this script. The script allows you to use different noise levels,
% different types of noise (uniforma and gaussian) and to systematically
% change localization error in case of CMD algoirthms.
%
% for a simpler test script that vary a single dimensions see testcmd1chg.m
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, G-SteX: Greedy Stem Extension for
% Free-Length Constrained Motif Discovery, IEA/AIE 2012
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%
doFinalizations=0;

noises=[0,0.05,0.1];       
        % this is the noise to be added to the signal itself. 
        
locNoises=[-1];
        %    this is the noise to be added to the localization 
        %   (a fraction of  the motif length). This is only effective if
        %   usersst==0
        %   if you use a negative number rsst will be used for this
        %   locNoise value. For example locNoises=[0,.1,-1] will calculate
        %   using localization noises of 0 and .1 and then using rsst


usersst=1;
generateNew=1;
generateNewNoise=0;
generateNewLocNoise=0;
generateNewCPDLocs=0;
generateNewDGRResults=1;
applyCMD=1;
applyMDFInalization=0;
generateNewStats=0;
generateNewFinalizationStats=0;
trySlidingDuringFinalExtension=0;
showResults=1;
delayList=0:5:25;
nSignals=2;
nMotifsToUse=2;
nOccursMax=6;
nOccursMin=3;
minMotif=60; 
maxMotif=60; 
M=10; 
K=12;
strLen=minMotif;
nSymbols=10;
nAlphabit=20;
nIterations=ceil(.7*nSymbols);
fractionKept=0.7;

algNames={'GEMODA+SAX';'GEMODA Euclidean'; 'Constrained GEMODA' ;'Shift Density'};
algFuns={@gmoda;@gemoda;@gmoda;@sdCMD};
algUsesCP=[0,0     ,1    ,1];
algParams={...
        {[]};
        {[]};
        {[]};{[]}...
        };

if ~doFinalizations
    generateNewFinalizationStats=0;
    applyMDFInalization=0;
end
autoParams=0;
firstW=-0.5;
threshold=0.05;
maxpcon=inf;
useBalancedStats=1;
minOverlapFraction=0;


if usersst==1
    locNoises=0;
end

calcrsst=0;
% we can add rsst to the list of noises we calculate
if min(locNoises)<0
    calcrsst=1;
end

nAlgs=size(algFuns,1);
nNoises=length(noises);
nLocNoises=length(locNoises);
nPoints=length(delayList);
locs=cell(nSignals,nNoises,nLocNoises,nAlgs);
means=cell(nSignals,nNoises,nLocNoises,nAlgs);
stats=cell(nSignals,nNoises,nLocNoises,nAlgs);
times=zeros(nSignals,nNoises,nLocNoises,nAlgs);

locsF=cell(nSignals,nNoises,nLocNoises,nAlgs);
meansF=cell(nSignals,nNoises,nLocNoises,nAlgs);
statsF=cell(nSignals,nNoises,nLocNoises,nAlgs);
timesF=zeros(nSignals,nNoises,nLocNoises,nAlgs);

locsT=cell(nSignals,1);
meansT=cell(nSignals,1);
statsT=cell(nSignals,1);

signals=cell(nSignals,nNoises+1);
strings=cell(nSignals,nNoises+1);
strTimes=cell(nSignals,nNoises+1);
cpdTimes=zeros(nSignals,nNoises+1);

if generateNew
    generateNewNoise=1;    
    generateNewLocNoise=1;    
end
if generateNewNoise || generateNewLocNoise
    generateNewCPDLocs=1;    
end
if generateNewCPDLocs
    generateNewDGRResults=1;    
end
if generateNewDGRResults
    generateNewStats=1;
end

save('./Data/params.mat');

if generateNew      
    for s=1:nSignals        
        [x,motifsT{s},locsT{s},occursT{s}]=generatePatterns(...
            nMotifsToUse,randi([nOccursMin,nOccursMax],1,nMotifsToUse)...
            ,[maxMotif,maxMotif],1,[-0.5,0.5]...
            ,4*maxMotif,0,0);    
        signals{s,1}=x;        
        tttt=tic();        
        strings{s,1}=timeseries2symbol(x,strLen,nSymbols,nAlphabit);
        strTimes{s,1}=toc(tttt);
    end    
    save('./Data/signals.mat','signals','motifsT','locsT','occursT','M'...
        ,'K','minMotif','maxMotif','strLen','strings','strTimes',...
        'nSymbols','nIterations','fractionKept','doFinalizations');
else
    load('./Data/params.mat');
    load('./Data/signals.mat','signals','motifsT','locsT','occursT'...
        ,'M','K','minMotif','maxMotif','strLen','strings','strTimes','nSymbols'...
        ,'nIterations','fractionKept','doFinalizations');
    nSignals=size(signals,1);    
    save('./Data/params.mat');
end

if generateNewDGRResults
    for s=1:nSignals
        % remove motifs with single occurrences if any (redundant!!)
        n=length(locsT{s});
        for i=1:n
            no=size(locsT{s}{i},1);
            if no<2
                motifsT(i)=[];
                locsT{s}(i)=[];
                occursT(i)=[];
            end
        end

        if(autoParams)
            n=length(locsT{s});
            if(n>0)
                mn=inf; mx=-inf;
                for i=1:n
                    mn2=min(locsT{s}{i}(:,2)-locsT{s}{i}(:,1)+1);
                    mx2=min(locsT{s}{i}(:,2)-locsT{s}{i}(:,1)+1);
                    if mn2<mn; mn=mn2; end;
                    if mx2>mx; mx=mx2; end;
                end
                minMotif=max(2,floor(mn*.9));
                maxMotif=max(2,ceil(mx*1.1));
                M=max(2,floor(minMotif/2));
                K=round(M*1.1);
            end
        end

        x=signals{s,1};
        xmax=max(x); xmin=min(x);
        T=numel(x);
        noiseA=cell(nNoises,1);
        if generateNewNoise
            for n=1:nNoises            
                noiseA{n}=noises(n).*0.5.*(xmax-xmin).*randn(size(x));                              
            end
            save(sprintf('./Data/noise_%d.mat',s),'noiseA');
        else
            load(sprintf('./Data/noise_%d.mat',s),'noiseA');
        end   
        for n=1:nNoises
            y=x+noiseA{n};
            signals{s,n+1}=y;
            
            tttt=tic();
            yD=timeseries2symbol(y,strLen,nSymbols,nAlphabit);
            strings{s,n+1}=yD;
            strTimes{s,n+1}=toc(tttt);
            if usersst
                if generateNewCPDLocs
                    tttt=tic();
                    [pcon]=cpd(y,M,K,-1,'sst3defaults',1,'dfunname',...
                        'distAvgEigDist','normalize',1,'calchankelevery',0,...
                        'npts4h',0,'symnorm',0,'scaleOutput',1);
                    cpdTimes(s,n)=toc(tttt);
                    pcon=pcon./max(pcon); 
                    pcon(pcon>maxpcon)=maxpcon; 
                    pcon=pcon./max(pcon);
                    candLocs=findLocsTh(pcon,threshold,firstW);
                    candLocs(candLocs>max(size(x))-minMotif+1)=[];
                    save(sprintf('./Data/rsstLocs_%03d.mat',s),...
                        'candLocs');                
                else
                    load(sprintf('./Data/rsstLocs_%03d.mat',s),...
                        'candLocs');                
                end
            end
            for l=1:nLocNoises
            %[locs{s,n,a},occurs]=cmd(x,20,105,0,pcon);
                if applyCMD
                     a=1; strt=tic;
                    fprintf('%s(Signal %d of %d, Noise %d of %d, Loc Noise %d of %d, Alg.  %d of %d) with motif range [%d,%d]...'...
                                    ,algNames{1},s,nSignals,n,nNoises,l,nLocNoises,1,nAlgs...
                                    ,strLen,...
                                    strLen);
                    [locs{s,n,l,1}] = gemoda(yD,...
                                        minMotif,2,[],@hamming,[],[],[],[],maxMotif);
                                    fprintf('DONE\n');
                    if ~algUsesCP(a)
                            if l==1
                                times(s,n,l,a)=toc(strt)/T+strTimes{s,n+1};
                            else
                                times(s,n,l,a)=times(s,n,l-1,a);
                            end
                        else
                            times(s,n,l,a)=toc(strt)/T+cpdTimes(s,n);
                    end
                       a=2; strt=tic;
                    
                    fprintf('%s(Signal %d of %d, Noise %d of %d, Loc Noise %d of %d, Alg.  %d of %d) with motif range [%d,%d]...'...
                                    ,algNames{2},s,nSignals,n,nNoises,l,nLocNoises,2,nAlgs...
                                    ,strLen,...
                                    strLen);
                    [locs{s,n,l,2}] = gemoda(y,...
                                        minMotif,0.1,[],[],[],[],[],[],maxMotif);
                                    fprintf('DONE\n');
                    if ~algUsesCP(a)
                            if l==1
                                times(s,n,l,a)=toc(strt)/T+strTimes{s,n+1};
                            else
                                times(s,n,l,a)=times(s,n,l-1,a);
                            end
                        else
                            times(s,n,l,a)=toc(strt)/T+cpdTimes(s,n);
                    end
                      a=3;  strt=tic;
                    
                   fprintf('%s(Signal %d of %d, Noise %d of %d, Loc Noise %d of %d, Alg.  %d of %d) with motif range [%d,%d]...'...
                                    ,algNames{3},s,nSignals,n,nNoises,l,nLocNoises,3,nAlgs...
                                    ,strLen,...
                                    strLen);
                                
                    [locs{s,n,l,3}] = gemoda(y,...
                                        minMotif,0.1,[],[],[],[],candLocs,[],maxMotif);
                                    fprintf('DONE\n');
                    if ~algUsesCP(a)
                            if l==1
                                times(s,n,l,a)=toc(strt)/T+strTimes{s,n+1};
                            else
                                times(s,n,l,a)=times(s,n,l-1,a);
                            end
                        else
                            times(s,n,l,a)=toc(strt)/T+cpdTimes(s,n);
                    end
                    
                else
                    load(sprintf('./Data/results.mat'),'times','locs','means','stats');
                end
                for a=4:nAlgs                           
                    if applyCMD
                        if ~algUsesCP(a)
                            if l==1
                                fprintf('%s(Signal %d of %d, Noise %d of %d, Loc Noise %d of %d, Alg.  %d of %d) with motif range [%d,%d]...'...
                                ,algNames{a},s,nSignals,n,nNoises,l,nLocNoises,a,nAlgs...
                                ,strLen,...
                                strLen); 
                                strt=tic;
                            else
                                fprintf('%s(Signal %d of %d, Noise %d of %d, Loc Noise Do not care, Alg.  %d of %d) with motif range [%d,%d]...'...
                                ,algNames{a},s,nSignals,n,nNoises,a,nAlgs...
                                ,strLen,...
                                strLen);                                 
                            end
                        else
                            fprintf('%s(Signal %d of %d, Noise %d of %d, Loc Noise %d of %d, Alg.  %d of %d) with motif range [%d,%d]...'...
                            ,algNames{a},s,nSignals,n,nNoises,l,nLocNoises,a,nAlgs...
                            ,minMotif,...
                            maxMotif); 
                            strt=tic;
                        end
                        if ~algUsesCP(a)
                            if l==1
                                [locs{s,n,l,a},means{s,n,l,a},stats{s,n,l,a}] = algFuns{a}(y,...
                                [minMotif,maxMotif]);
                            else
                                locs{s,n,l,a}=locs{s,n,l-1,a};
                                means{s,n,l,a}=means{s,n,l-1,a};
                                stats{s,n,l,a}=stats{s,n,l-1,a};
                            end
                        else                            
                            [locs{s,n,l,a},means{s,n,l,a},stats{s,n,l,a}] = algFuns{a}(y,...
                                [minMotif,maxMotif],candLocs);
                        end
                        if ~algUsesCP(a)
                            if l==1
                                times(s,n,l,a)=toc(strt)/T+strTimes{s,n+1};
                            else
                                times(s,n,l,a)=times(s,n,l-1,a);
                            end
                        else
                            times(s,n,l,a)=toc(strt)/T+cpdTimes(s,n);
                        end
                        fprintf(' Done\n');
                    else
                        load(sprintf('./Data/results.mat'),'times','locs','means','stats');
                    end
                    if applyMDFInalization
                        fprintf('Finalization(Signal %d of %d, Noise %d of %d, Loc Noise %d of %d, Alg.  %d of %d) with motif range [%d,%d]...'...
                        ,s,nSignals,n,nNoises,l,nLocNoises,a,nAlgs...
                        ,minMotif,...
                        maxMotif); 
                        strt=tic;
                        [locsF{s,n,l,a},meansF{s,n,l,a},statsF{s,n,l,a}] = mdFinalize(y,[minMotif,maxMotif],...
                            candLocs,locs{s,n,l,a},means{s,n,l,a},stats{s,n,l,a},'trySlidingDuringFinalExtension',trySlidingDuringFinalExtension);
                        timesF(s,n,l,a)=toc(strt)/T;
                        fprintf(' Done\n');
                    end
                    
                end 
            end
        end
    end
    if applyCMD
        save(sprintf('./Data/results.mat'),'times','locs','means','stats');
        save('./Data/signals.mat','signals','motifsT','locsT','occursT','M'...
        ,'K','minMotif','maxMotif','strLen','strings','strTimes','nSymbols','nIterations','fractionKept');
    else
        load(sprintf('./Data/results.mat'),'times','locs','means','stats');        
    end
    if applyMDFInalization
        save(sprintf('./Data/resultsF.mat'),'timesF','locsF','meansF','statsF');
    else
        load(sprintf('./Data/results.mat'),'timesF','locsF','meansF','statsF');
    end
else
    load(sprintf('./Data/results.mat'),'times','locs','means','stats');
    if doFinalizations
        load(sprintf('./Data/resultsF.mat'),'timesF','locsF','meansF','statsF');
    end
end


if generateNewStats
    statsCPD=cell(nSignals,nNoises,nLocNoises,nAlgs,nPoints);
    delayStats=cell(nSignals,nNoises,nLocNoises,nAlgs,nPoints);
    statsA=cell(nSignals,nNoises,nLocNoises,nAlgs);
    statsM=cell(nSignals,nNoises,nLocNoises,nAlgs);        
    coveringNone=zeros(nSignals,nNoises,nLocNoises,nAlgs);
    coveringPartial=zeros(nSignals,nNoises,nLocNoises,nAlgs);
    coveringMulti=zeros(nSignals,nNoises,nLocNoises,nAlgs);
    covered=zeros(nSignals,nNoises,nLocNoises,nAlgs,nMotifsToUse);
    extras=zeros(nSignals,nNoises,nLocNoises,nAlgs,nMotifsToUse);
    badMotif=cell(nSignals,nNoises,nLocNoises,nAlgs);
    covering=cell(nSignals,nNoises,nLocNoises,nAlgs);
    coveredby=cell(nSignals,nNoises,nLocNoises,nAlgs);
    for s=1:nSignals
        T=numel(signals{s});
        for n=1:nNoises
            for l=1:nLocNoises
                for a=1:nAlgs
                    fprintf('Stats(Signal %d of %d, Noise %d of %d, Loc Noise %d of %d, Alg.  %d of %d) with motif range [%d,%d]...'...
                        ,s,nSignals,n,nNoises,l,nLocNoises,a,nAlgs...
                        ,minMotif,...
                        maxMotif); 
                    [foundLocs,trueLocs]=mdq2cpq(locs{s,n,l,a},locsT{s});
                    statsA{s,n,l,a}=mdq(locs{s,n,l,a},locsT{s},T);
                    statsM{s,n,l,a}=mdqM(locs{s,n,l,a},locsT{s},T,minOverlapFraction);
                    for d=1:nPoints
                        [statsCPD{s,n,l,a,d},delayStats{s,n,l,a,d}]=...
                            cpquality(foundLocs,trueLocs...
                            ,T,nPoints,useBalancedStats,@findLocsTh);
                    end
                    
                    [coveringNone(s,n,l,a),...
                        coveringPartial(s,n,l,a),...
                        coveringMulti(s,n,l,a),...
                        covered(s,n,l,a,:),...
                        extras(s,n,l,a,:),...
                        badMotif{s,n,l,a},covering{s,n,l,a},coveredby{s,n,l,a}]...
                        =mdquality(locs{s,n,l,a},locsT{s});
                    fprintf(' Done\n');
                end
            end
        end
    end
    statsCPD=cell2mat(statsCPD);
    statsA=cell2mat(statsA);
    statsM=cell2mat(statsM);
    delayStats=cell2mat(delayStats);
    save('./Data/stats.mat','statsA','statsM','statsCPD','delayStats'...
        ,'coveringNone','coveringPartial','coveringMulti','covered','extras','badMotif'...
        ,'covering','coveredby');
else
    load('./Data/stats.mat','statsA','statsM','statsCPD','delayStats'...
        ,'coveringNone','coveringPartial','coveringMulti','covered','extras','badMotif'...
        ,'covering','coveredby');
end

if generateNewFinalizationStats    && doFinalizations
    statsCPDF=cell(nSignals,nNoises,nLocNoises,nAlgs,nPoints);
    delayStatsF=cell(nSignals,nNoises,nLocNoises,nAlgs,nPoints);
    statsAF=cell(nSignals,nNoises,nLocNoises,nAlgs);
    statsMF=cell(nSignals,nNoises,nLocNoises,nAlgs);        
    coveringNoneF=zeros(nSignals,nNoises,nLocNoises,nAlgs);
    coveringPartialF=zeros(nSignals,nNoises,nLocNoises,nAlgs);
    coveringMultiF=zeros(nSignals,nNoises,nLocNoises,nAlgs);
    coveredF=zeros(nSignals,nNoises,nLocNoises,nAlgs,nMotifsToUse);
    extrasF=zeros(nSignals,nNoises,nLocNoises,nAlgs,nMotifsToUse);
    badMotifF=cell(nSignals,nNoises,nLocNoises,nAlgs);
    coveringF=cell(nSignals,nNoises,nLocNoises,nAlgs);
    coveredbyF=cell(nSignals,nNoises,nLocNoises,nAlgs);
    for s=1:nSignals
        T=numel(signals{s});
        for n=1:nNoises
            for l=1:nLocNoises
                for a=1:nAlgs
                    fprintf('Final Stats(Signal %d of %d, Noise %d of %d, Loc Noise %d of %d, Alg.  %d of %d) with motif range [%d,%d]...'...
                        ,s,nSignals,n,nNoises,l,nLocNoises,a,nAlgs...
                        ,minMotif,...
                        maxMotif); 
                    [foundLocsF,trueLocsF]=mdq2cpq(locsF{s,n,l,a},locsT{s});
                    statsAF{s,n,l,a}=mdq(locsF{s,n,l,a},locsT{s},T);
                    statsMF{s,n,l,a}=mdqM(locsF{s,n,l,a},locsT{s},T,minOverlapFraction);
                    for d=1:nPoints
                        [statsCPDF{s,n,l,a,d},delayStatsF{s,n,l,a,d}]=...
                            cpquality(foundLocsF,trueLocsF...
                            ,T,nPoints,useBalancedStats,@findLocsTh);
                    end
                    
                    [coveringNoneF(s,n,l,a),...
                        coveringPartialF(s,n,l,a),...
                        coveringMultiF(s,n,l,a),...
                        coveredF(s,n,l,a,:),...
                        extrasF(s,n,l,a,:),...
                        badMotifF{s,n,l,a},coveringF{s,n,l,a},coveredbyF{s,n,l,a}]...
                        =mdquality(locsF{s,n,l,a},locsT{s});
                    fprintf(' Done\n');
                end
            end
        end
    end
    
    statsCPDF=cell2mat(statsCPDF);
    statsAF=cell2mat(statsAF);
    statsMF=cell2mat(statsMF);
    delayStatsF=cell2mat(delayStatsF);
    save('./Data/statsF.mat','statsAF','statsMF','statsCPDF','delayStatsF'...
        ,'coveringNoneF','coveringPartialF','coveringMultiF','coveredF','extrasF','badMotifF'...
        ,'coveringF','coveredbyF');    
else
    if (doFinalizations)
    load('./Data/statsF.mat','statsAF','statsMF','statsCPDF','delayStatsF'...
        ,'coveringNoneF','coveringPartialF','coveringMultiF','coveredF','extrasF','badMotifF'...
        ,'coveringF','coveredbyF');
    end
end
if showResults
    displayTestResults;
end
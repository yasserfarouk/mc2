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

% Tests differnt constraint motif discovery algorithms
% See testcmdLNoises.m for the parameters that need to be set before using
% this script
%
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
addpath('../cpd');
%addpath('../Projections');
%addpath('./sax');
addpath('../common/matlab');
doFinalizations=0;

useLocNoiseAsLength=1;
nLengths=5;

generateNew=1;
generateNewNoise=1;
generateNewLocNoise=1;
generateNewCPDLocs=1;
generateNewDGRResults=1;
applyCMD=1;
applyMDFInalization=doFinalizations;
generateNewStats=1;
generateNewFinalizationStats=1;
trySlidingDuringFinalExtension=0;
showResults=1;
delayList=0:5:25;
nSignals=5;
nMotifsToUse=4;
nOccursMax=4;
nOccursMin=2;
nMotifLengthMaxBetweenMotifs=3;
minMotif=25; 
maxMotif=75;
lStepNormal=min((maxMotif-minMotif),max(4,ceil((maxMotif-minMotif)/10)));
noises=[0];%,0.001,0.01,0.1];       
        % this is the noise to be added to the signal itself. 
        
locNoises=[0];%[0,0.000001,0.00001,0.001,0.1];

usersst=0;

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
algNames={'NMK+','MK++','MK++ Recalc','MK++ Recalc All'}; %,'MK--','MK-- Recalc','MK-- Recalc All'
algFuns={@mkpRange;@mkppRange;@mkppRange;@mkppRange};%;@mkmmRange;@mkmmRange;@mkmmRange
nMotifs=10.*[10,10,10,10,10,10,10];
lStep=.5*[10,10,10,10,10,10,10];
finalParam=[0,0,1,2,0,1,2]; % will be recalcDist in case of mkpp and Fast in case of mkp


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
            ,nMotifLengthMaxBetweenMotifs*maxMotif,0,0);    
        signals{s,1}=x;        
        tttt=tic();        
        %strings{s,1}=timeseries2symbol(x,strLen,nSymbols,nAlphabit);
        %strTimes{s,1}=toc(tttt);
    end    
    save('./Data/signals.mat','signals','motifsT','locsT','occursT'...
        ,'minMotif','maxMotif',...
        'doFinalizations');
else
    load('./Data/params.mat');
    load('./Data/signals.mat','signals','motifsT','locsT','occursT'...
        ,'minMotif','maxMotif','doFinalizations');
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
        sPerm=randperm(nSignals);
        for l1=1:nLengths
            if useLocNoiseAsLength
                x=[];
                nX=0;
                for cc=1:l1
                    x=[x,signals{sPerm(l1),1}];                    
                end
            else
                x=signals{s,1};                
            end            
            xmax=max(x); xmin=min(x);
            T=numel(x);
            noiseA=cell(nNoises,1);
            if generateNewNoise
                for n=1:nNoises            
                    noiseA{n}=noises(n).*0.5.*(xmax-xmin).*randn(size(x));                              
                end
                if l1==1
                    save(sprintf('./Data/noise_%d.mat',s),'noiseA');
                else
                    save(sprintf('./Data/noise_%d_%d.mat',s,l1),'noiseA');
                end
            else
                if l1==1
                    load(sprintf('./Data/noise_%d.mat',s),'noiseA');
                else
                    load(sprintf('./Data/noise_%d_%d.mat',s,l1),'noiseA');
                end
            end 
            for n=1:nNoises
                y=x+noiseA{n};
                if l1==1
                    signals{s,n+1}=y;
                end
                for l2=1:nLocNoises
                %[locs{s,n,a},occurs]=cmd(x,20,105,0,pcon);
                    l=max(l1,l2);
                    for a=1:nAlgs                    
                        if applyCMD
                            fprintf('%s(Signal %d of %d, Noise %d of %d, Loc Noise %d of %d, Alg.  %d of %d) with motif range [%d,%d]...'...
                                    ,algNames{a},s,nSignals,n,nNoises,l,nLocNoises,a,nAlgs...
                                    ,minMotif,...
                                    maxMotif); 
                                    strt=tic;
                                    try
                                        timeIStarted=tic;
                                        [locs{s,n,l,a},means{s,n,l,a},stats{s,n,l,a},dummy,...
                                            times(s,n,l,a)] = ...
                                            algFuns{a}(y,[minMotif,maxMotif],lStep(a),...
                                            nMotifs(a),finalParam(a));                        
                                        times(s,n,l,a)=times(s,n,l,a)/T;
                                        fprintf(sprintf(' Done in %d us/point\n',times(s,n,l,a)));
                                    catch e                                    
                                        times(s,n,l,a)=inf;
                                        fprintf(sprintf(' FAILED in %d us\n',1e6*toc(timeIStarted)));
                                    end
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
                    if useLocNoiseAsLength
                        break;
                    end
                end 
            end
        end
    end
    if applyCMD
        save(sprintf('./Data/results.mat'),'times','locs','means','stats');
        save('./Data/signals.mat','signals','motifsT','locsT','occursT'...
            ,'minMotif','maxMotif');
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
    load('./Data/statsF.mat','statsAF','statsMF','statsCPDF','delayStatsF'...
        ,'coveringNoneF','coveringPartialF','coveringMultiF','coveredF','extrasF','badMotifF'...
        ,'coveringF','coveredbyF');
end

if showResults
    displayTestResults;
end
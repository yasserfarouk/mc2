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
% this script. The script allows you to use different noise levels,
% different types of noise (uniforma and gaussian) and to systematically
% change motif lengths.
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
if ~doFinalizations
    generateNewFinalizationStats=0;
    applyMDFInalization=0;
end
if ~exist('minMotif2Gen','var')
    minMotif2Gen=minMotif;
end
if ~exist('maxMotif2Gen','var')
    maxMotif2Gen=maxMotif;
end
autoParams=0;
firstW=-0.5;
threshold=0.05;
maxpcon=inf;
useBalancedStats=1;
minOverlapFraction=0;

nAlgs=size(algFuns,1);
nNoises=length(noises);
nLengths=length(lengths);
nPoints=length(delayList);
locs=cell(nSignals,nNoises,nLengths,nAlgs);
means=cell(nSignals,nNoises,nLengths,nAlgs);
stats=cell(nSignals,nNoises,nLengths,nAlgs);
times=nan(nSignals,nNoises,nLengths,nAlgs);

locsF=cell(nSignals,nNoises,nLengths,nAlgs);
meansF=cell(nSignals,nNoises,nLengths,nAlgs);
statsF=cell(nSignals,nNoises,nLengths,nAlgs);
timesF=nan(nSignals,nNoises,nLengths,nAlgs);

locsT=cell(nSignals,1);
meansT=cell(nSignals,1);
statsT=cell(nSignals,1);

signals=cell(nSignals,nNoises+1,nLengths);
strings=cell(nSignals,nNoises+1,nLengths);
%strTimes=cell(nSignals,nNoises+1,nLengths);
cpdTimes=zeros(nSignals,nNoises+1,nLengths);

if generateNew
    generateNewNoise=1;        
end
if generateNewNoise 
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
            ,[minMotif2Gen,maxMotif2Gen],1,[-0.5,0.5]...
            ,nMotifsBetween*maxMotif2Gen,0,0);    
        signals{s,1,1}=x;        
%         tttt=tic();        
%         strings{s,1,1}=timeseries2symbol(x,strLen,nSymbols,nAlphabit);
%         strTimes{s,1,1}=toc(tttt);
    end
    nMax=0;
    for s=1:nSignals
        lo=cell2mat(locsT{s});
        lo=lo(:,2);
        nMax=max([lo;nMax]);
    end
    nMax=nMax+maxMotif2Gen;
    lengths=lengths.*nMax;
    save('./Data/signals.mat','signals','motifsT','locsT','occursT','M'...
        ,'K','minMotif','maxMotif','strLen','strings',...
        'nSymbols','nIterations','fractionKept','doFinalizations');
else
    load('./Data/params.mat');
    load('./Data/signals.mat','signals','motifsT','locsT','occursT'...
        ,'M','K','minMotif','maxMotif','strLen','strings','nSymbols'...
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

        x=signals{s,1,1};
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
            for l=1:nLengths                
                [y,locsT{s},~]=embedInside(y,locsT{s},y,lengths(l),[-0.5,0.5],0,false,true);
                signals{s,n+1,l}=y;
                if generateNewCPDLocs
                    if max(algUsesCP)>0.5
                        tttt=tic();
                        [pcon]=cpd(y,M,K,-1,'sst3defaults',1,'dfunname',...
                            'distAvgEigDist','normalize',1,'calchankelevery',0,...
                            'npts4h',0,'symnorm',0,'scaleOutput',1);
                        cpdTimes(s,n,l)=toc(tttt);
                        pcon=pcon./max(pcon); 
                        pcon(pcon>maxpcon)=maxpcon; 
                        pcon=pcon./max(pcon);
                        candLocs=findLocsTh(pcon,threshold,firstW);
                        save(sprintf('./Data/rsstLocs_%03d_%03d_%03d.mat',s,n,l),...
                            'candLocs');     
                    end
                else
                    load(sprintf('./Data/rsstLocs_%03d_%03d_%03d.mat',s,n,l),...
                        'candLocs');                
                end                                       
                                                  
                for a=1:nAlgs                                                    
                    if applyCMD
                        if 0
                            fprintf('%s(Signal %d of %d, Noise %d of %d, Length %d of %d=%d, Alg.  %d of %d) with motif range [%d,%d]...'...
                                ,algNames{a},s,nSignals,n,nNoises,l,nLengths,lengths(l),a,nAlgs...
                                ,strLen,...
                                strLen); 
                                strt=tic;
                        else
                            fprintf('%s(Signal %d of %d, Noise %d of %d, Length %d of %d=%d, Alg.  %d of %d) with motif range [%d,%d]...'...
                            ,algNames{a},s,nSignals,n,nNoises,l,nLengths,lengths(l),a,nAlgs...
                            ,minMotif,...
                            maxMotif); 
                            strt=tic;
                        end
%                         try
                            strt=tic;
                            if ~isempty(algParams{a}{1})
                                if ~algUsesCP(a)
                                    if strcmp(tolower(algNames{a}),'projections')
                                        yD=timeseries2symbol(y,strLen,nSymbols,nAlphabit);
                                        strings{s,n+1,l}=yD;
                                        %strTimes{s,n+1,l}=toc(tttt);  
                                        [locs{s,n,l,a},means{s,n,l,a},stats{s,n,l,a}] = algFuns{a}(yD,...
                                                nIterations,fractionKept,y,strLen/nSymbols,[],algParams{a});
                                    else
                                        [locs{s,n,l,a},means{s,n,l,a},stats{s,n,l,a}] = algFuns{a}(y,[minMotif,maxMotif]);                                        
                                    end
                                else
                                    [locs{s,n,l,a},means{s,n,l,a},stats{s,n,l,a}] = algFuns{a}(y,...
                                        [minMotif,maxMotif],candLocs,algParams{a});
                                end
                            else
                                 if ~algUsesCP(a)
                                    if strcmpi((algNames{a}),'projections')
                                        [locs{s,n,l,a},means{s,n,l,a},stats{s,n,l,a}] = algFuns{a}(yD,...
                                            nIterations,fractionKept,y,strLen/nSymbols);
                                    else
                                        [locs{s,n,l,a},means{s,n,l,a},stats{s,n,l,a}] = algFuns{a}(y,[minMotif,maxMotif]);
                                    end
                                else
                                    [locs{s,n,l,a},means{s,n,l,a},stats{s,n,l,a}] = algFuns{a}(y,...
                                        [minMotif,maxMotif],candLocs);
                                 end                                                        
                            end                            
                            times(s,n,l,a)=toc(strt);
                            if algUsesCP(a)
                                times(s,n,l,a)=times(s,n,l,a)+cpdTimes(s,n,l);
                            end
                            times(s,n,l,a)=times(s,n,l,a)/T;
                            fprintf(sprintf(' Done(%d motifs found)\n',numel(locs{s,n,l,a})));
%                         catch e
%                             fprintf(' FAILED\n');
%                         end
                        
                    else
                        load(sprintf('./Data/results.mat'),'times','locs','means','stats');
                        
                    end
                    if applyMDFInalization
                        fprintf('Finalization(Signal %d of %d, Noise %d of %d, Length %d of %d=%d, Alg.  %d of %d) with motif range [%d,%d]...'...
                        ,s,nSignals,n,nNoises,l,nLengths,lengths(l),a,nAlgs...
                        ,minMotif,...
                        maxMotif); 
                        strt=tic;
                        [locsF{s,n,l,a},meansF{s,n,l,a},statsF{s,n,l,a}] = mdFinalize(y,[minMotif,maxMotif],...
                            candLocs,locs{s,n,l,a},means{s,n,l,a},stats{s,n,l,a},'trySlidingDuringFinalExtension',trySlidingDuringFinalExtension);
                        timesF(s,n,l,a)=toc(strt)/T;
                        fprintf(sprintf(' Done(%d motifs found)\n',numel(locsF{s,n,l,a})));
                    end
                    
                end 
                
            end
        end
    end
    if applyCMD
        save(sprintf('./Data/results.mat'),'times','locs','means','stats');
        save('./Data/signals.mat','signals','motifsT','locsT','occursT','M'...
        ,'K','minMotif','maxMotif','strLen','strings','nSymbols','nIterations','fractionKept');        
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
    statsCPD=cell(nSignals,nNoises,nLengths,nAlgs,nPoints);
    delayStats=cell(nSignals,nNoises,nLengths,nAlgs,nPoints);
    statsA=cell(nSignals,nNoises,nLengths,nAlgs);
    statsM=cell(nSignals,nNoises,nLengths,nAlgs);        
    coveringNone=zeros(nSignals,nNoises,nLengths,nAlgs);
    coveringPartial=zeros(nSignals,nNoises,nLengths,nAlgs);
    coveringMulti=zeros(nSignals,nNoises,nLengths,nAlgs);
    covered=zeros(nSignals,nNoises,nLengths,nAlgs,nMotifsToUse);
    extras=zeros(nSignals,nNoises,nLengths,nAlgs,nMotifsToUse);
    badMotif=cell(nSignals,nNoises,nLengths,nAlgs);
    covering=cell(nSignals,nNoises,nLengths,nAlgs);
    coveredby=cell(nSignals,nNoises,nLengths,nAlgs);
    for s=1:nSignals        
        for n=1:nNoises
            for l=1:nLengths
                %T=numel(signals{s,n,l});
                T=lengths(l);
                for a=1:nAlgs
                    fprintf('Stats(Signal %d of %d, Noise %d of %d, Loc Noise %d of %d, Alg.  %d of %d) with motif range [%d,%d]...'...
                        ,s,nSignals,n,nNoises,l,nLengths,a,nAlgs...
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
    statsCPDF=cell(nSignals,nNoises,nLengths,nAlgs,nPoints);
    delayStatsF=cell(nSignals,nNoises,nLengths,nAlgs,nPoints);
    statsAF=cell(nSignals,nNoises,nLengths,nAlgs);
    statsMF=cell(nSignals,nNoises,nLengths,nAlgs);        
    coveringNoneF=zeros(nSignals,nNoises,nLengths,nAlgs);
    coveringPartialF=zeros(nSignals,nNoises,nLengths,nAlgs);
    coveringMultiF=zeros(nSignals,nNoises,nLengths,nAlgs);
    coveredF=zeros(nSignals,nNoises,nLengths,nAlgs,nMotifsToUse);
    extrasF=zeros(nSignals,nNoises,nLengths,nAlgs,nMotifsToUse);
    badMotifF=cell(nSignals,nNoises,nLengths,nAlgs);
    coveringF=cell(nSignals,nNoises,nLengths,nAlgs);
    coveredbyF=cell(nSignals,nNoises,nLengths,nAlgs);
    for s=1:nSignals
        
        for n=1:nNoises
            for l=1:nLengths
                %T=numel(signals{s,n,l});
                T=lengths(l);
                for a=1:nAlgs
                    fprintf('Final Stats(Signal %d of %d, Noise %d of %d, Loc Noise %d of %d, Alg.  %d of %d) with motif range [%d,%d]...'...
                        ,s,nSignals,n,nNoises,l,nLengths,a,nAlgs...
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
    if(doFinalizations)
    load('./Data/statsF.mat','statsAF','statsMF','statsCPDF','delayStatsF'...
        ,'coveringNoneF','coveringPartialF','coveringMultiF','coveredF','extrasF','badMotifF'...
        ,'coveringF','coveredbyF');
    end
end
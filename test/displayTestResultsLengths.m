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

% displays the results of testcmd
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%
load('./Data/params.mat');
load('./Data/signals.mat');
load('./Data/results.mat');
load('./Data/stats.mat');
if doFinalizations
load('./Data/resultsF.mat');
load('./Data/statsF.mat');
end



showSignalsAndMotifs=1;
useBeforeFinalization=1-doFinalizations;
showFinalizationResults=doFinalizations;
showRSST=0;
noSubplots=1;
showCPD=max(algUsesCP);

if ~showFinalizationResults
    useBeforeFinalization=1;
end
signalsToDisp=1:nSignals;
noisesToDisp=1:nNoises;
lengthsToDisp=1;%1:nLengths;
algsToDisp=1:nAlgs;


colors={'r','b','g','k','c','m','--r','--b','--g','--k','--c','--m'...
,'-.r','-.b','-.g','-.k','-.c','-.m',':r',':b',':g',':k',':c',':m'};

nAlgsToDisp=numel(algsToDisp);
if showSignalsAndMotifs
    for s=signalsToDisp        
        %x=signals{s}; T=numel(x);
        locsTCurrent=locsT{s};
        motifsTCurrent=motifsT{s};
        %load(sprintf('./Data/lnoise_%03d.mat',s));
        load(sprintf('./Data/noise_%d.mat',s));
        for n=noisesToDisp
            %y=x+noiseA{n};
            for l=lengthsToDisp 
                x=signals{s,1+1,l};
                y=signals{s,n+1,l};
                
                T=numel(y);
               figure;
               subplot(nAlgsToDisp+1,1,nAlgsToDisp+1);                           
                
                plot(y,':y');  
                if showCPD
                    load(sprintf('./Data/rsstLocs_%03d.mat',s),...
                            'candLocs'); 
                    locv=zeros(size(y));
                    locv(candLocs)=max(y);
                    locv(locv<0.5)=inf;
                    hold on; plot(locv,'*r'); 
                end    
                
                hold off;
                nMotifs=length(locsTCurrent);
                hold on;
                for i=1:nMotifs    
                    tmp=zeros(size(y)); tmp(tmp<1)=inf;
                    clr=colors{mod(i,length(colors))+1};
                    for j=1:size(locsTCurrent{i},1);
                        tmp(locsTCurrent{i}(j,1):locsTCurrent{i}(j,2))=...
                            y(locsTCurrent{i}(j,1):locsTCurrent{i}(j,2));
                    end
                    plot(tmp,clr);
                end                    
                title(sprintf('True Motif Locations(Signal:%d)',s));
                if showRSST
                    %subplot(nRows,2,[nxt+2:nxt+3]);
                    subplot(nAlgsToDisp+1,1,nAlgsToDisp+1);
                    hold on;
                    plot(pcon+ones(size(pcon))); 
                end
                for a=algsToDisp
                    if useBeforeFinalization
                        meansCurrent=means{s,n,l,a};
                        locsCurrent=locs{s,n,l,a};
                        statsCurrent=stats{s,n,l,a};
                    else
                        meansCurrent=meansF{s,n,l,a};
                        locsCurrent=locsF{s,n,l,a};
                        statsCurrent=statsF{s,n,l,a};
                    end
                    mx=max(y)+.05.*max(y); mn=min(y)-.05.*min(y);

                    subplot(nAlgsToDisp+1,1,a); 
                    plot(y,':y');  
                    if showCPD
                        locv=zeros(size(y));
                        locv(candLocs)=max(y);
                        locv(locv<0.5)=inf;
                        hold on; plot(locv,'*r'); 
                    end
                    hold off;
                    nMotifs=length(meansCurrent);
                    colors={'r','b','g','k','c','m','--r','--b','--g','--k','--c','--m'...
                        ,'-.r','-.b','-.g','-.k','-.c','-.m',':r',':b',':g',':k',':c',':m'};
                    hold on;
                    nMotifs=numel(locsCurrent);
                    for i=1:nMotifs    
                        tmp=zeros(size(y)); tmp(tmp<1)=inf;
                        clr=colors{mod(i,length(colors))+1};
                        for j=1:size(locsCurrent{i},1);
                            tmp(locsCurrent{i}(j,1):locsCurrent{i}(j,2))=y(locsCurrent{i}(j,1):locsCurrent{i}(j,2));
                        end
                        plot(tmp,clr);
                    end
                    
                    title(sprintf('%s(Signal:%d,Noise:%f,Loc Noise:%f)',algNames{a},s,noises(n),lengths( l)));                                       
                end
                
            end
        end
    end
end


coveringTrue=1-coveringMulti-coveringNone-coveringPartial;
if showFinalizationResults
coveringTrueF=1-coveringMultiF-coveringNoneF-coveringPartialF;
end

figure1=figure; 
covPerAlg=[];
for a=1:nAlgs
    %if a==3; continue; end;
    covPerAlgt=[];
    tmp=coveringTrue(:,:,1,a);
    covPerAlgt=[covPerAlgt,mean(reshape(tmp,numel(tmp),1))];    
    tmp=coveringPartial(:,:,1,a);
    covPerAlgt=[covPerAlgt,mean(reshape(tmp,numel(tmp),1))];
    tmp=coveringMulti(:,:,1,a);
    covPerAlgt=[covPerAlgt,mean(reshape(tmp,numel(tmp),1))];
    tmp=coveringNone(:,:,1,a);
    covPerAlgt=[covPerAlgt,mean(reshape(tmp,numel(tmp),1))];
    covPerAlg=[covPerAlg;covPerAlgt];
end
%figure
h=subplot(1+showFinalizationResults,2,1);
bar(covPerAlg,'BarLayout','stacked','BarWidth',0.2);
set(h,'XTickLabel',algNames,...
    'XTick',[1 2 3 4 5]);

legend({'Correct','Covering Partial', 'Covering Multiple','Covering None'});
title('Results Before Finalization');

covPerAlg=[];
for a=1:nAlgs
    %if a==3; continue; end;
    covPerAlgt=[];
    tmp=covered(:,:,1,a,:);
    covPerAlgt=[covPerAlgt,mean(reshape(tmp,numel(tmp),1))];
    tmp=extras(:,:,1,a,:);
    covPerAlgt=[covPerAlgt,mean(reshape(tmp,numel(tmp),1))];
    covPerAlg=[covPerAlg;covPerAlgt];    
end
%figure
h=subplot(1+showFinalizationResults,2,2);
bar(covPerAlg,'BarLayout','grouped');
set(h,'XTickLabel',algNames,...
    'XTick',[1 2 3 4 5]);

legend({'Covered','Extras'});
title('Results Before Finalization');


%figure; 
if showFinalizationResults
covPerAlg=[];
for a=1:nAlgs
   %if a==3; continue; end;
    covPerAlgt=[];
    tmp=coveringTrueF(:,:,1,a);
    covPerAlgt=[covPerAlgt,mean(reshape(tmp,numel(tmp),1))];    
    tmp=coveringPartialF(:,:,1,a);
    covPerAlgt=[covPerAlgt,mean(reshape(tmp,numel(tmp),1))];
    tmp=coveringMultiF(:,:,1,a);
    covPerAlgt=[covPerAlgt,mean(reshape(tmp,numel(tmp),1))];
    tmp=coveringNoneF(:,:,1,a);
    covPerAlgt=[covPerAlgt,mean(reshape(tmp,numel(tmp),1))];
    covPerAlg=[covPerAlg;covPerAlgt];
end
%figure
h=subplot(2,2,3);
bar(covPerAlg,'BarLayout','stacked','BarWidth',0.2);
set(h,'XTickLabel',algNames,...
    'XTick',[1 2 3 4 5]);

legend({'Correct','Covering Partial', 'Covering Multiple','Covering None'});
title('Results After Finalization');

covPerAlg=[];
for a=1:nAlgs
    %if a==3; continue; end;
    covPerAlgt=[];
    tmp=coveredF(:,:,1,a,:);
    covPerAlgt=[covPerAlgt,mean(reshape(tmp,numel(tmp),1))];
    tmp=extrasF(:,:,1,a,:);
    covPerAlgt=[covPerAlgt,mean(reshape(tmp,numel(tmp),1))];
    covPerAlg=[covPerAlg;covPerAlgt];    
end
%figure
h=subplot(2,2,4);
bar(covPerAlg,'BarLayout','grouped');
set(h,'XTickLabel',algNames,...
    'XTick',[1 2 3 4 5]);

legend({'Covered','Extras'});
title('Results After Finalization');
end

% display the effect of noise on the true coverage
figure;
hold on;
for a=1:nAlgs
    tmp=squeeze(coveringTrue(:,1,:,a));
    mymean=mean(tmp,1);
    mystd=std(tmp,1);            
    plot(lengths,mymean,colors{a},'LineWidth',2);
    set(gca,'XScale','log','FontSize',24);
end        
legend(algNames);
ylabel('Correct Discovery Rate');


% display the effect of localization noise on the covered
figure;
hold on;
for a=1:nAlgs
    tmp=squeeze(mean(squeeze((covered(:,1,:,a,:))),3));
    mymean=mean(tmp,1);
    mystd=std(tmp,1);            
    plot(lengths,mymean,colors{a},'LineWidth',2);
    set(gca,'XScale','log','FontSize',24);
end        
legend(algNames);
ylabel('Fraction of True Motifs Covered');

% display the effect of noise on the extra
figure;
hold on;
for a=1:nAlgs
    tmp=squeeze(mean(squeeze((extras(:,1,:,a,:))),3));
    mymean=mean(tmp,1);
    mystd=std(tmp,1);            
    plot(lengths,mymean,colors{a},'LineWidth',2);
    set(gca,'XScale','log','FontSize',24);
end        
legend(algNames);
ylabel('Fraction of Extra Parts in Discovered Motifs');


%     tmp=squeeze(mean(squeeze((covered(:,:,1,a,:))),3));
%     mymean=mean(tmp,1);
%     mystd=std(tmp,1);    
%     if noSubplots
%        figure(h2);               
%     else
%         h=subplot(1,3,2,'XScale','log');
%     end
%     hold on;
%     plot(noises,mymean,colors{a},'LineWidth',2);
%     
%     tmp=squeeze(mean(squeeze((extras(:,:,1,a,:))),3));
%     mymean=mean(tmp,1);
%     mystd=std(tmp,1);    
%     if noSubplots
%        figure(h3);        
%     else
%         h=subplot(1,3,3,'XScale','log');
%     end
%     hold on;
%     plot(noises,mymean,colors{a},'LineWidth',2);
% end
% statNames={'Covering True','Covered','Extra'};
% if noSubplots
%     hh=[h1,h2,h3];
%     for i=1:3
%         figure(hh(i));
%         a1=axes(hh(i));
%         set(a1,'XScale','log');
%         legend(algNames);
%         title(sprintf('Efect of noise on %s value',statNames{i}));    
%     end
% else
%     for i=1:3    
%         h=subplot(1,3,3,'XScale','log');
%         legend(algNames);
%         title(sprintf('Efect of noise on %s value',statNames{i}));    
%     end
% end

if showFinalizationResults
    % display the effect of noise on the true coverage
    figure;
    hold on;
    for a=1:nAlgs
        tmp=squeeze(coveringTrueF(:,1,:,a));
        mymean=mean(tmp,1);
        mystd=std(tmp,1);            
        plot(lengths,mymean,colors{a},'LineWidth',2);
        set(gca,'XScale','log','FontSize',24);
    end        
    legend(algNames);
    ylabel('Correct Discovery Rate');


    % display the effect of localization noise on the covered
    figure;
    hold on;
    for a=1:nAlgs
        tmp=squeeze(mean(squeeze((coveredF(:,1,:,a,:))),3));
        mymean=mean(tmp,1);
        mystd=std(tmp,1);            
        plot(lengths,mymean,colors{a},'LineWidth',2);
        set(gca,'XScale','log','FontSize',24);
    end        
    legend(algNames);
    ylabel('Fraction of True Motifs Covered');

    % display the effect of noise on the extra
    figure;
    hold on;
    for a=1:nAlgs
        tmp=squeeze(mean(squeeze((extrasF(:,1,:,a,:))),3));
        mymean=mean(tmp,1);
        mystd=std(tmp,1);            
        plot(lengths,mymean,colors{a},'LineWidth',2);
        set(gca,'XScale','log','FontSize',24);
    end        
    legend(algNames);
    ylabel('Fraction of Extra Parts in Discovered Motifs');
end

allTimes=zeros(size(times,1)*size(times,2)*size(times,3),size(times,4));
if showFinalizationResults
allTimesF=zeros(size(times,1)*size(times,2)*size(times,3),size(times,4));
end
for a=1:nAlgs
    nxt=1;
    for s=1:nSignals
        for n=1:nNoises
            for l=1:nLengths
                allTimes(nxt,a)=times(s,n,l,a);
                if showFinalizationResults
                    allTimesF(nxt,a)=timesF(s,n,l,a);
                end
                nxt=nxt+1;
            end
        end
    end
end
figure;
boxplot(allTimes,'notch','on');
set(gca,'XTickLabel',algNames,...
    'XTick',[1 2 3 4 5]);
ylabel('Time Without Finalization');

if showFinalizationResults
figure;
boxplot(allTimesF+allTimes,'notch','on');
set(gca,'XTickLabel',algNames,...
    'XTick',[1 2 3 4 5]);

ylabel('Time with Finalization');
end
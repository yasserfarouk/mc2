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

% display the results of running testcpd
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, On Comparing SSA-based Change Point
% Discovery Algorithms, IEEE SII 2011, 938-945 
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%
clinespec={'--r','--g','--b','--c','--m','--y','--k','r','g','b','c','m','y','k',':r',':g',':b',':c',':m',':y',':k',':r',':g',':b',':c',':m',':y',':k'};
clinespecFew={'r','--g',':b','*c','om',':k'};
displayControlFig=0;
displayESRN=0;
displayESRL=0;
nMaxPerPage=9;    
statsDisp=5:8;%1:nStats
chToDisplay=1:nAlgs;%[6,7,35,15,23,31];    
nPoints4Curves=1;

    close all;
    figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for j=1:numel(statsDisp)
        i=statsDisp(j);
        subplot(ceil(numel(statsDisp)/2),2,j);
        hold on;
        for current=1:nAlgs
            tmp=squeeze(statVals(:,nTh+1,:,i,current));
            if(showErrorBars)
                errorbar(points,mean(tmp,1),std(tmp,1),clinespec{mod(current,28)+1});             
            else
                plot(points,mean(tmp,1),clinespec{mod(current,28)+1});             
            end
        end
        title(sprintf('mean %s',statNames{i}));
    end
    %legend(chNames); 
    if saveFigs
        savefig('./Figs/stats',gcf, 'png','eps');
    end
if displayControlFig    
    figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    for j=1:numel(statsDisp)
        i=statsDisp(j);
        subplot(ceil(numel(statsDisp)/2)+1,2,j);
        hold on;
        for current=1:6
            tmp=squeeze(statVals(:,nTh+1,:,i,current));
            if(showErrorBars)
                errorbar(points,mean(tmp,1),std(tmp,1),clinespecFew{mod(current,length(clinespecFew))+1});             
            else
                plot(points,mean(tmp,1),clinespecFew{mod(current,length(clinespecFew))+1});             
            end
        end
        title(sprintf('mean %s',statNames{i}));
    end
    subplot(ceil(numel(statsDisp)/2)+1,2,numel(statsDisp)+1);
    hold on;
    for current=1:6
        tmp=squeeze(esrs(:,current,:));
        if(showErrorBars)
            errorbar(points,mean(tmp,1),std(tmp,1),clinespecFew{mod(current,length(clinespecFew))+1});
        else
            plot(points,mean(tmp,1),clinespecFew{mod(current,length(clinespecFew))+1});             
        end
    end
    title(sprintf('mean ESR'));
if displayESRN    
    subplot(ceil(numel(statsDisp)/2)+1,2,numel(statsDisp)+2);
    hold on;
    for current=1:6
        tmp=squeeze(esrse(:,current,1:nP));
        if(showErrorBars)
            errorbar(points,mean(tmp,1),std(tmp,1),clinespecFew{mod(current,length(clinespecFew))+1});
        else
            plot(points,mean(tmp,1),clinespecFew{mod(current,length(clinespecFew))+1});             
        end
    end
    title(sprintf('mean ESR(Exponential Scaling)'));
end
    legend(chNames(1:6),'Location','NorthWestOutside');
    %legend(chNames); 
    if saveFigs
        savefig('./Figs/controlalgs',gcf, 'png','eps');
    end
end    
    figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    subplot(1,2,1);
    hold on;    
    for current=1:nAlgs        
        tmp=squeeze(esrs(:,current,:));
        if(showErrorBars)
            errorbar(points,mean(tmp,1),std(tmp,1),clinespec{mod(current,28)+1});             
        else
            plot(points,mean(tmp,1),clinespec{mod(current,28)+1});             
        end
    end
    title('ESR');
if displayESRN    
    subplot(1,2,2);
    hold on;    
    for current=2:nAlgs        
        tmp=squeeze(esrse(:,current,1:nP));
        if(showErrorBars)
            errorbar(points,mean(tmp,1),std(tmp,1),clinespec{mod(current,28)+1});             
        else
            plot(points,mean(tmp,1),clinespec{mod(current,28)+1});             
        end
    end
    title('ESR (Exponential Scaling)');    
end
    if saveFigs
        savefig('./Figs/esr2',gcf, 'png','eps');
    end
    
    
    figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,1);
    hold on;
    for current=1:nAlgs
        tmp=squeeze(cpsps(:,current,:));
        if(showErrorBars)
            errorbar(points,mean(tmp,1),std(tmp,1),clinespec{mod(current,28)+1});             
        else
            plot(points,mean(tmp,1),clinespec{mod(current,28)+1});             
        end
    end
    title('CPSP');
    legend(chNames,'Location','NorthWestOutside');    
    %legend(chNames,'Location','NorthWest'); 
    hold off;

    subplot(2,2,2);
    hold on;    
    for current=1:nAlgs        
        tmp=squeeze(esrs(:,current,:));
        if(showErrorBars)
            errorbar(points,mean(tmp,1),std(tmp,1),clinespec{mod(current,28)+1});             
        else
            plot(points,mean(tmp,1),clinespec{mod(current,28)+1});             
        end
    end
    title('ESR');
    %legend(chNames); 
    

    subplot(2,2,3);
   
    hold on;
    for current=2:nAlgs    
        tmp=squeeze(desrs(:,current,:));
        if(showErrorBars)
            errorbar(points,mean(tmp,1),std(tmp,1),clinespec{mod(current,28)+1});             
        else
            plot(points,mean(tmp,1),clinespec{mod(current,28)+1});             
        end
    end
    title('ESR Difference');
    

    subplot(2,2,4);
    hold on;
    for current=1:nAlgs                
        tmp=squeeze(aesrs(:,current,:));
        if(showErrorBars)
            errorbar(points,mean(tmp,1),std(tmp,1),clinespec{mod(current,28)+1});             
        else
            plot(points,mean(tmp,1),clinespec{mod(current,28)+1});             
        end
    end
    title('Accomulated ESR');
    %legend(chNames); 
    if saveFigs
    savefig('./Figs/esr',gcf, 'png','eps');
    end
    
if displayESRL    
    figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,1);
    hold on;    
    for current=1:nAlgs        
        tmp=squeeze(esrs(:,current,:));
        if(showErrorBars)
            errorbar(points,mean(tmp,1),std(tmp,1),clinespec{mod(current,28)+1});             
        else
            plot(points,mean(tmp,1),clinespec{mod(current,28)+1});             
        end
    end
    title('ESR');
    legend(chNames,'Location','NorthWestOutside');    
    %legend(chNames,'Location','NorthWest'); 
    hold off;

    subplot(2,2,2);
    hold on;    
    for current=1:nAlgs        
        tmp=squeeze(esrsl(:,current,:,nTh+1));
        if(showErrorBars)
            errorbar(points,mean(tmp,1),std(tmp,1),clinespec{mod(current,28)+1});             
        else
            plot(points,mean(tmp,1),clinespec{mod(current,28)+1});             
        end
    end
    title('ESR (After Localization)');
    %legend(chNames); 


    subplot(2,2,3);
   
    hold on;
    for current=2:nAlgs    
        tmp=squeeze(desrsl(:,current,:));
        if(showErrorBars)
            errorbar(points,mean(tmp,1),std(tmp,1),clinespec{mod(current,28)+1});             
        else
            plot(points,mean(tmp,1),clinespec{mod(current,28)+1});             
        end
    end
    title('ESR Difference (After Localization)');
    

    subplot(2,2,4);
    hold on;
    for current=1:nAlgs                
        tmp=squeeze(aesrsl(:,current,:));
        if(showErrorBars)
            errorbar(points,mean(tmp,1),std(tmp,1),clinespec{mod(current,28)+1});             
        else
            plot(points,mean(tmp,1),clinespec{mod(current,28)+1});             
        end
    end
    title('Accomulated ESR (After Localization)');
    %legend(chNames); 
    if saveFigs
    savefig('./Figs/esrl',gcf, 'png','eps');
    end
    
    figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    subplot(1,2,1);
    hold on;    
    for current=1:nAlgs        
        tmp=squeeze(esrsl(:,current,:,nTh+1));
        if(showErrorBars)
            errorbar(points,mean(tmp,1),std(tmp,1),clinespec{mod(current,28)+1});             
        else
            plot(points,mean(tmp,1),clinespec{mod(current,28)+1});             
        end
    end
    title('ESR (After Localization)');
        
    bxs=[];
    sp=sum(1./[1:nP]);
    for current=1:nAlgs                
        tmp=squeeze(esrsl(:,current,:,nTh+1));
        for s=1:nSeries
            for p=1:nP
                tmp(s,p)=tmp(s,p)./(p);
            end
        end
        bxs=[bxs,sum(tmp,2)];
    end
    bxs=bxs./sp;
    %axes1 = axes('Parent',h,'XTickLabel',chNames);    
    axes1=subplot(1,2,2);
    %box(axes1,'on');
    %hold(axes1,'all');
    boxplot(bxs,'plotstyle','compact');
    set(axes1,'XTick',1:nAlgs);
    %set(axes1,'XTickLabel',chNames(:)');
    title('AUAESR (After Localization)');
    
    if saveFigs
        savefig('./Figs/esr3',gcf, 'png','eps');
    end
end
if displayESRN    
    figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    axes1=subplot(2,2,1);
    set(axes1,'XTick',1:nAlgs);
    set(axes1,'XTickLabel',chNames);
    box(axes1,'on');
    hold(axes1,'all');
    tmp=squeeze(esrse(:,:,nP+1));
    boxplot(tmp,'plotstyle','compact');    
    title('ESR (Exponential Scaling) with INF allowable delay');    
    %


    subplot(2,2,2);
    hold on;    
    for current=2:nAlgs        
        tmp=squeeze(esrse(:,current,1:nP));
        if(showErrorBars)
            errorbar(points,mean(tmp,1),std(tmp,1),clinespec{mod(current,28)+1});             
        else
            plot(points,mean(tmp,1),clinespec{mod(current,28)+1});             
        end
    end
    title('ESR (Exponential Scaling)');    
    

    subplot(2,2,3);  
    hold on;
    for current=2:nAlgs    
        tmp=squeeze(desrse(:,current,1:nP));
        if(showErrorBars)
            errorbar(points,mean(tmp,1),std(tmp,1),clinespec{mod(current,28)+1});             
        else
            plot(points,mean(tmp,1),clinespec{mod(current,28)+1});             
        end
    end
    title('ESR Difference (Exponential Scaling)');
    legend(chNames(2:nAlgs),'Location','SouthWestOutside'); 

    subplot(2,2,4);
    hold on;
    for current=2:nAlgs                
        tmp=squeeze(aesrse(:,current,1:nP));
        if(showErrorBars)
            errorbar(points,mean(tmp,1),std(tmp,1),clinespec{mod(current,28)+1});             
        else
            plot(points,mean(tmp,1),clinespec{mod(current,28)+1});             
        end
    end
    title('Accomulated ESR (Exponential Scaling)');
    %legend(chNames); 
    if saveFigs
    savefig('./Figs/esre',gcf, 'png','eps');
    end
end  %savex(x{series});
    
    figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);    
    bxs=[];
    sp=sum(1./[1:nP]);
    for current=1:nAlgs                
        tmp=squeeze(cpsps(:,current,:));
        for s=1:nSeries
            for p=1:nP
                tmp(s,p)=tmp(s,p)./(p);
            end
        end
        bxs=[bxs,sum(tmp,2)];
    end
    bxs=bxs./sp;
    %axes1 = axes('Parent',h,'XTickLabel',chNames);        
if displayESRN
    axes1=subplot(2,1,1);    
    %box(axes1,'on');
    %hold(axes1,'all');
    tmp=squeeze(esrse(:,:,nP+1));
    boxplot(tmp,'plotstyle','compact');    
    set(axes1,'XTick',1:nAlgs);
    %set(axes1,'XTickLabel',chNames(:)');
    title('ESR (Exponential Scaling) with INF allowable delay');    
%     set(axes1,'XTick',1:nAlgs);
%     set(axes1,'XTickLabel',chNames);
%     box(axes1,'on');
%     hold(axes1,'all');
%     boxplot(bxs,'plotstyle','compact');
%     title('Discrete Version of ESR');
end       
    bxs=[];
    sp=sum(1./[1:nP]);
    for current=1:nAlgs                
        tmp=squeeze(esrs(:,current,:));
        for s=1:nSeries
            for p=1:nP
                tmp(s,p)=tmp(s,p)./(p);
            end
        end
        bxs=[bxs,sum(tmp,2)];
    end
    bxs=bxs./sp;
    %axes1 = axes('Parent',h,'XTickLabel',chNames);    
    axes1=subplot(2,1,2);
    %box(axes1,'on');
    %hold(axes1,'all');
    boxplot(bxs,'plotstyle','compact');
    set(axes1,'XTick',1:nAlgs);
    %set(axes1,'XTickLabel',chNames(:)');
    title('ESR Score using ESR(p)/p');
    
%     h=subplot(2,2,3);
%     bxs=[];
%     for current=1:nAlgs                
%         tmp=squeeze(desrs(:,current,:));        
%         for s=1:nSeries
%             for p=1:nP
%                 tmp(s,p)=tmp(s,p)./(p);
%             end            
%         end        
%         bxs=[bxs,sum(tmp,2)];
%     end    
%     %axes1 = axes('Parent',h,'XTickLabel',chNames);    
%     axes1=h;
%     set(axes1,'XTick',1:nAlgs);
%     set(axes1,'XTickLabel',chNames);
%     box(axes1,'on');
%     hold(axes1,'all');
%     boxplot(bxs,'plotstyle','compact');
%     title('ESR Score using DESR(p)/p');
        
%     h=subplot(2,2,3);
%     bxs=[];
%     for current=1:nAlgs                
%         tmp=squeeze(aesrs(:,current,:));
%         bxs=[bxs,sum(tmp,2)];
%     end
%     bxs=bxs./nP;
%     %axes1 = axes('Parent',h,'XTickLabel',chNames);    
%     axes1=h;
%     set(axes1,'XTick',1:nAlgs);
%     set(axes1,'XTickLabel',chNames);
%     box(axes1,'on');
%     hold(axes1,'all');
%     boxplot(bxs,'plotstyle','compact');
%     title('ESR Score using sum(AESR(p))/nP');
%     
    h=figure;
    axes1 = axes('Parent',h,'XTickLabel',chNames);    
    set(axes1,'XTick',1:nAlgs);
    set(axes1,'XTickLabel',chNames(1:end));
    box(axes1,'on');
    hold(axes1,'all');
    boxplot(1000*times,'plotstyle','compact');
    title('Execution Time (ms/point)');
    if saveFigs
    savefig('./Figs/time',gcf, 'png','eps');
    end
% evaluate performance curves ROC etc    
    
    %statVals=zeros(nSeries,nTh+1,nP,nStats,nAlgs);    
    figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,1);
    hold on;
    sp=sum(1./(1:nP));
    for current=1:nAlgs
        x1=1-squeeze(statVals(:,:,:,find(ismember(statNames,'Specificity')==1),current));
        y1=squeeze(statVals(:,:,:,find(ismember(statNames,'Sensitivity')==1),current));
        x1=squeeze(mean(x1,1)); y1=squeeze(mean(y1,1));
        for p=1:nP
            x1(:,p)=x1(:,p)./p;
            y1(:,p)=y1(:,p)./p;
        end        
        x1=squeeze(sum(x1,2)); y1=squeeze(sum(y1,2));
        x1=x1./sp; y1=y1./sp;
        [x1,x1O]=sort(x1); y1=y1(x1O);
        plot(x1, y1,clinespec{mod(current,28)+1}); 
    end
    axis([0,1,0,1]);
    title('ROC (sensitivity vs. (1-specificity))');
    legend(chNames,'Location','NorthWestOutside');
    %legend(chNames); 

    subplot(2,2,2);
    hold on;
    for current=1:nAlgs
        x1=squeeze(statVals(:,:,:,find(ismember(statNames,'Precision')==1),current));
        y1=squeeze(statVals(:,:,:,find(ismember(statNames,'Recall')==1),current));
        x1=squeeze(mean(x1,1)); y1=squeeze(mean(y1,1));
        for p=1:nP
            x1(:,p)=x1(:,p)./p;
            y1(:,p)=y1(:,p)./p;
        end        
        x1=squeeze(sum(x1,2)); y1=squeeze(sum(y1,2));
        x1=x1./sp; y1=y1./sp;
        [x1,x1O]=sort(x1); y1=y1(x1O);
        plot(x1, y1, clinespec{mod(current,28)+1}); 
    end
    axis([0,1,0,1]);
    title('recall vs. prec');
    %legend(chNames); 

    subplot(2,2,3);
    hold on;
    for current=1:nAlgs
        x1=squeeze(statVals(:,:,:,find(ismember(statNames,'Accuracy')==1),current));
        y1=squeeze(statVals(:,:,:,find(ismember(statNames,'Specificity')==1),current));
        x1=squeeze(mean(x1,1)); y1=squeeze(mean(y1,1));
        for p=1:nP
            x1(:,p)=x1(:,p)./p;
            y1(:,p)=y1(:,p)./p;
        end        
        x1=squeeze(sum(x1,2)); y1=squeeze(sum(y1,2));
        x1=x1./sp; y1=y1./sp;
        [x1,x1O]=sort(x1); y1=y1(x1O);
        plot(x1, y1, clinespec{mod(current,28)+1}); 
    end
    axis([0,1,0,1]);
    title('specificity vs. accuracy');    
    if saveFigs
    savefig('./Figs/curves',gcf, 'png','eps');
    end
    
% display best results (latest series)    
    if( displaySeriesNum>0)    
        series=displaySeriesNum;
        ch=changes{series};
        chLocs=changeLocs{series};
        nBegin=1;
        tmpt=t{series}; tmpt(tmpt<0.0001)=nan;
        while(1)
            figure;
            set(gcf,'units','normalized','outerposition',[0 0 1 1]);
            todisp=chToDisplay(nBegin:min(nBegin+nMaxPerPage,length(chToDisplay))); 
            nDisp=length(todisp)+1;
            subplot(nDisp,1,1); plot(x{series}); hold on; 
            plot(t{series},'r'); hold off;
            %subplot(nAlgs+2,1,2); %plot(xTrue{series}); hold on; plot(t{series},'r'); hold off;
            cplot=2;            
            for current=todisp
                subplot(nDisp,1,cplot);  cplot=cplot+1;
                plot(tmpt,'or'); hold on; 
                plot(ch(current,:),'b'); 

%                 nBest=thBest(series,nPoints4Curves, ...
%                     find(ismember(statNames,'MCC')==1),current);
%                 tmp=zeros(size(t{series})); tmp(chLocs{current,nBest})=1; tmp(tmp==0)=nan;
%                 plot(tmp,'or'); 

%                 nBestF=thBest(series,nPoints4Curves, ...
%                     find(ismember(statNames,'F')==1),current);
%                 tmp=zeros(size(t{series})); tmp(chLocs{current,nBestF})=0.8; tmp(tmp==0)=nan;
%                 plot(tmp,'*g'); 
                %ylabel(sprintf('%s\nTh(MCC)=%g\nTh(F1)=%g',chNames{current,1},thresholds(nBest),thresholds(nBestF))); 
%                 [mxesr,nBestE]=max(squeeze(esrsl(series,current,:,1:nTh)),[],2);                
%                 tmp=zeros(size(t{series})); tmp(chLocs{current,nBestE(nPoints4Curves)})=0.8; tmp(tmp==0)=nan;
%                 plot(tmp,'*g'); 
%                ylabel(sprintf('%s\nTh=%g\nESR=%g',chNames{current,1},thresholds(nBestE(nPoints4Curves)),mxesr(nPoints4Curves))); 
                ylabel(sprintf('%s',chNames{current,1})); 
                hold off; 
            end            
            if saveFigs
            savefig(sprintf('./Figs/outputs_%d',nBegin),gcf, 'png','eps');
            end
            if(nBegin+nMaxPerPage>=length(chToDisplay))
                break;
            end
        end
    end
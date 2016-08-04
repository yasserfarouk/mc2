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

close all;
clrs={'r';'b';'g';'m';'c';'y';'k'}; nclrs=length(clrs);
specs={'-';'--';':';'-.'};  nspecs=length(specs);
mrkrs={'+';'o';'x';'.';'*';'s';'d';'v';'^';'<';'>';'p';'h'}; nmrkrs=length(mrkrs);
detVis=clrs;   ndetvis=nclrs;             %change this to specs to make a single color plot

%save everything
save('./Data/syntheticAll.mat');

%display time delay stats.
figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
for meth=1:nm
    for det=1:nd
        if(~applicableDet(meth,det))
            continue;
        end
        
    end
end

% display results for every statistic (no averaging)
finalStatF=cell(nF,nm,nd);
figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
for fld=1:nF    
    subplot(nF,1,fld);
    hold on;
    lgd=cell(nd*nm,1);
    nxtLgd=1;    
    for meth=1:nm
        for det=1:nd        
            if(~applicableDet(meth,det))
                continue;
            end
            avgStatF=statF{meth,det,1,fld};
            for model=2:nModels        
                avgStatF=avgStatF+statF{meth,det,model,fld};
            end
            avgStatF=avgStatF./nModels;
            finalStatF{fld,meth,det}=avgStatF;
            plot(th,avgStatF,[specs{mod(meth-1,nspecs)+1},detVis{mod(det-1,ndetvis)+1}]);
            lgd{nxtLgd}=sprintf('%s (%s)',method{meth},detector{det});
            nxtLgd=nxtLgd+1;
        end
    end
    toremove=[];
    for l=1:length(lgd)
        if(isempty(lgd{l}))
            toremove=[toremove,l];
        end
    end
    lgd(toremove)=[];
    legend(lgd,'Location','SouthEast');
    title (fields(fld));    
    hold off;
end
savefig('./Figs/stats',gcf,'Parent', 'png','eps');
% display box plot of every statistic
figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
for fld=1:nF
    subplot(nF,1,fld);
    lgd=cell(nd*nm,1);
    nxtLgd=1;
    boxs=[];
    for meth=1:nm
        for det=1:nd        
            if(~applicableDet(meth,det))
                continue;
            end
            boxs=[boxs,finalStatF{fld,meth,det}(:)];            
            lgd{nxtLgd}=sprintf('%s (%s)',method{meth},detector{det});
            nxtLgd=nxtLgd+1;
        end
    end    
    toremove=[];
    for l=1:length(lgd)
        if(isempty(lgd{l}))
            toremove=[toremove,l];
        end
    end
    lgd(toremove)=[];
    boxplot(boxs);
    set(gca,'XTickLabel',lgd,'XTick',1:length(lgd));
    title (fields(fld));    
end
savefig('./Figs/statsbox',gcf, 'png','eps');

% display box plot of times
figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
lgd=cell(nd*nm,1);
nxtLgd=1;
boxs=[];
for meth=1:nm
    for det=1:nd        
        if(~applicableDet(meth,det))
            continue;
        end
        if(nModels==1)
            boxs=[boxs,[mytimes(meth,det,:);mytimes(meth,det,:)]];            
        else
            boxs=[boxs,squeeze(mytimes(meth,det,:))];            
        end
        lgd{nxtLgd}=sprintf('%s (%s)',method{meth},detector{det});
        nxtLgd=nxtLgd+1;
    end
end    
toremove=[];
for l=1:length(lgd)
    if(isempty(lgd{l}))
        toremove=[toremove,l];
    end
end
lgd(toremove)=[];
boxplot(boxs);
set(gca,'XTickLabel',lgd,'XTick',1:length(lgd));
title ('Execution Time Per Point'); 
savefig('./Figs/timebox',gcf, 'png','eps');


% roc curve
figure;
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
lgd=cell(nd*nm,1);
nxtLgd=1;
hold on;
for meth=1:nm
    for det=1:nd                    
        if(~applicableDet(meth,det))
            continue;
        end
        mrcg1=rocg{meth,det,1}(:,1);
        mrcg2=rocg{meth,det,1}(:,2);
        for i=2:nModels
            mrcg1=mrcg1+rocg{meth,det,i}(:,1);
            mrcg2=mrcg2+rocg{meth,det,i}(:,2);
        end
        mrcg1=mrcg1./nModels; mrcg2=mrcg2./nModels;
        plot(mrcg2,mrcg1, ...
        [specs{mod(meth-1,nspecs)+1},detVis{mod(det-1,ndetvis)+1}]);                        
        lgd{nxtLgd}=sprintf('%s (%s)',method{meth},detector{det});
        nxtLgd=nxtLgd+1;
    end
end
axis([0,1,0,1]);
title('ROC curve');
toremove=[];
for l=1:length(lgd)
    if(isempty(lgd{l}))
        toremove=[toremove,l];
    end
end
lgd(toremove)=[];
legend(lgd,'Location','SouthEast');
hold off;
savefig('./Figs/roc',gcf, 'png','eps');

if(nModels<3)
    %display the signals and their changes
    for model=1:nModels        
        figure;
        set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        nx=size(c{model,det},1);
        for i=1:nx
            subplot(nx,2,2*(i-1)+1);
            plot(xs{model}(i,:));
            subplot(nx,2,2*(i));
            plot(chgss{model}(i,:),'r');    
            hold;
            for det=1:nd
                if~isempty(c{model,det}) && ~strcmpi(detector(det),'none') && ~strcmpi(detector(det),'optimal')
                      plot(c{model,det}(i,:)./max(c{model,det}(i,:)),[detVis{mod(det-1,ndetvis)+2}]);    
                end
            end
        end
    end
end
savefig('./Figs/signals',gcf, 'png','eps');
%display best structure according to every statistic for first model
for fld=1:nF
    figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    bestv=[];
    model=0;
    for modd=1:nModels            
        for meth=1:nm
            for det=1:nd
                if ~isempty(bestStat{meth,det,modd,fld}) && (isempty(bestv) || bestStat{meth,det,modd,fld}>bestv)
                    bestv=bestStat{meth,det,modd,fld};
                    model=modd;
                end
            end
        end
    end
    nProcesses=size (Ms{model},1);        
    subplot(nm+1,nd,nd+1);
    %drawModel(Ms{model},nProcesses);
    subplot(nm+1,nd,1);
    drawModel(Ms{model},nProcesses);
    title('Original Model');
%         if(sum(sum(abs(Es{model}-E2)))>0.00001)
%             warning('drawModel is not working correctly');
%         end
    
    for meth=1:nm
        for det=1:nd        
            if(~applicableDet(meth,det))
                continue;
            end
            subplot(nm+1,nd,(meth)*nd+det);            
            [dummy,Eg]=drawModel(MgSave{meth,det,model,fld},nProcesses);
            title(sprintf('%s(%s), threshold=%g, %s=%4.3f',method{meth} ...
                ,detector{det},thSaved{meth,det,model,fld},fields{fld}, ...
                bestStat{meth,det,model,fld}));
        end
    end    
    savefig(sprintf('./Figs/graphs%03d_%s',model,fields{fld}),gcf, 'png','eps');
end

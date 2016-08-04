function myboxplot(data,xlabels,ylbl,lineWidth,axesfontSize,ylblFontSize)
    if ~exist('ylblFontSize','var') || isempty(ylblFontSize)
        ylblFontSize=36;
    end        
    if ~exist('axesfontSize','var') || isempty(axesfontSize)
        axesfontSize=30;
    end
    if ~exist('lineWidth','var') || isempty(lineWidth)
        lineWidth=3;
    end 
    if(~exist('ylbl','var') || isempty(ylbl))
        ylbl='';
    end
    n=numel(xlabels);
    bh=boxplot(data,'notch','on');
    set(gca,'FontSize',axesfontSize);
    if exist('xlabels','var')
        set(bh,'linewidth',lineWidth);        
        set(gca,'XTick',1:n)
        set(gca,'XTickLabel',xlabels);    
    end
    if exist('ylbl','var')
        ylabel(ylbl,'FontSize',ylblFontSize);
    end
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
end
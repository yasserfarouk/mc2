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
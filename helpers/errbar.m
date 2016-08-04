function errbar(data,clr,err,xaxis,scale,xticks,xlbl,ylbl,lgnd,lineWidth,axesFontSize,xfontSize,yfontSize,lgndFontSize,lgndLoc)
% errbar(data,clr,err,xaxis,scale,xticks) plots error bars  
%
%  errbar(data,clr,err,xaxis,scale,xticks);  
%                               data - cell array of 2D data, err - 'std'/'stderr'  
%                               clr - colors vectors  
%                               err - 'std', 'stderr', 'mean'
%                               xaxis - optional [xmin:step:xmax] vector
%                               scale - 'log'  
%                               xticks - cellarray  
if ~exist('lineWidth','var') 
    lineWidth=3;
end

if ~exist('axesFontSize','var') 
    axesFontSize=30;
end
if ~exist('xfontSize','var') 
    xfontSize=axesFontSize+6;
end
if ~exist('yfontSize','var') 
    yfontSize=axesFontSize+6;
end
if exist('lgnd','var')
    if ~exist('lgndFontSize','var') 
        lgndFontSize=yfontSize;
    end
    if ~exist('lgndLoc','var') 
        lgndLoc='Best';
    end
end
if(nargin >= 4)
    x = xaxis;
    if(isempty(x))
        x = 1:size(data{1},2);
    end
    if(nargin >= 5)
        if(strcmpi(scale,'log'))
            for j=1:size(data,1)
                data{j} = 10*log10(data{j});
            end            
        end
    end    
else
    x = 1:size(data{1},2);
end


for j=1:size(data,1)
    if(length(x) > 1)
        dx(j) = 0.04 * (x(2) - x(1));        
    else
        dx(j) = 0.04;
    end
    if(mod(j,3) == 1)
        dx(j) = -dx(j);
    elseif(mod(j,3) == 0)
        dx(j) = 0;
    end
end

errBars=zeros(1,size(data,1));
for j=1:size(data,1)

    if(strcmpi(err,'std') == 1)
        e = std(data{j});
        errBars(j)=errorbar(x+dx(j),mean(data{j}),e,'-x','Color',clr(j,:));
        hold on;
    elseif(strcmpi(err,'stderr') == 1)
        e = std(data{j}) / size(data{j},1);
        errBars(j)=errorbar(x+dx(j),mean(data{j}),e,'-x','Color',clr(j,:));
        hold on;
    else
        errBars(j)=plot(x,mean(data{j}),'-x','Color',clr(j,:));
        hold on;
    end   
   
end

set(errBars,'LineWidth',lineWidth);
set(gcf,'Color',[1 1 1]);

if(length(x) <= 15)
    set(gca, 'XTickLabelMode','manual');
    set(gca, 'XTickMode','manual');
    set(gca, 'XTick', x);

    if(nargin >= 6)
        set(gca, 'XTickLabel', char(xticks));    
    else
        set(gca, 'XTickLabel', num2str(x'));    
    end    
end
set(gca,'FontSize',axesFontSize);
if exist('xlbl','var')
    xlabel(xlbl,'FontSize',xfontSize);
end
if exist('ylbl','var')
    ylabel(ylbl,'FontSize',yfontSize);
end
if exist('lgnd','var')
    legend(lgnd,'Location',lgndLoc,'FontSize',lgndFontSize);    
end

set(gcf,'units','normalized','outerposition',[0 0 1 1]);
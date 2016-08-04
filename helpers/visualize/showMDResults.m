function [figs ] = showMDResults(d, locs,means,stats,normalizeMeans, nRows,nCols,figTitle,maxMotifs )
%Displays a discovered motif and its occurrences
%   
showMean=0;
showLocs=1;
showStats=0;
showLocsInOriginal=1;
if ~exist('maxMotifs','var') || isempty(maxMotifs)
    maxMotifs=numel(locs);
else
    locs=locs(1:min([numel(locs),maxMotifs]));
    %means=means(1:maxMotifs);
    %stats=stats(1:maxMotifs);
end
if ~exist('normalizeMeans','var') || isempty(normalizeMeans)
    normalizeMeans=0;
end
showTitle=nargin>6;
if ~exist('showMean','var') || ~exist('means','var') || isempty(showMean) || isempty(means)
    showMean=0;
else 
    if isempty(means); showMean=0; end;
end
if ~exist('showStats','var') || ~exist('stats','var') || isempty(showStats) || isempty(stats)
    showStats=0;
else 
    if isempty(stats); showStats=0; end;
end
if ~exist('nRows','var') || isempty(nRows)
    nRows=min(5,length(locs));
    nCols=floor(min(nRows*2,length(locs))/nRows);    
elseif ~exist('nCols','var') || isempty(nCols)
    n=length(locs);
    nCols=min(2,ceil(n/nRows));
end
nRows=max(1,nRows);
nCols=max(1,nCols);
nPerPage=nRows*nCols;
nLocsTotal=size(locs,1);
nPages=ceil(nLocsTotal/nPerPage);
figs=cell(nPages,1);
T=length(d);
colors={'b:','g:','m:','c:','k:','b--','g--','m--','c--','k--','b-.','g-.','m-.','c-.','k-.'};
colorsSolid={'b','g','m','c','k','b:','g:','m:','c:','k:','b--','g--','m--','c--','k--','b-.','g-.','m-.','c-.','k-.'};
for p=1:nPages    
    i=(p-1)*nPerPage+1;
    nToGo=min(i+nPerPage,nLocsTotal+1)-i;
    figs{p}=figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);                
    if showTitle
        suptitle(sprintf('%s %d of %d',figTitle,p,nPages));
    end
    if showLocsInOriginal
        subplot(nRows+1,nCols,nRows*nCols+1:nRows*nCols+nCols);
        plot(d,'y:');
        hold on;
        leg=cell(nToGo+1,1);        
        leg{1}='timeseries';
        for j=1:nToGo; 
            current=i+j-1;
            leg{j+1}=sprintf('Motif %d',current);
            nLocs=size(locs{current},1);             
            tmp=inf*ones(size(d));                
            for k=1:nLocs;                 
                l1=locs{current}(k,1);
                l2=locs{current}(k,2);
                if l1>0 && l2<=T
                    tmp(l1:l2)=d(l1:l2);                    
                end
            end; 
            plot(tmp,colorsSolid{mod(j,length(colorsSolid))+1},'LineWidth',4);             
        end
        legend(leg);        
    end
    for j=1:nToGo; 
        current=i+j-1;
        subplot(nRows+showLocsInOriginal,nCols,j);  
        hold on; 
        nLocs=size(locs{current},1); 
        if showLocs
            for k=1:nLocs; 
                l1=locs{current}(k,1);
                l2=locs{current}(k,2);
                if l1>0 && l2<=T
                    if normalizeMeans
                        plot(normalize(d(l1:l2)),colors{mod(k,length(colors))+1}); 
                    else
                        plot(d(l1:l2),colors{mod(k,length(colors))+1}); 
                    end
                end
            end; 
        end
        if showMean
            plot(normalize(means{current}),'r','LineWidth',2); 
        end
        if showMean
            ylabel(sprintf('ARM %d \n|M| %d, |m| %d',...
                current,nLocs,length(means{current})));            
        else
            ylabel(sprintf('ARM %d \n|M| %d, |m| %d',...
                current,nLocs,(1+locs{current}(1,2)-locs{current}(1,1))));            
        end
        if showStats
            title(sprintf('mean %6.3f, var %6.3d',...
                stats(current).mean,stats(current).var));                    
        end
    end;
end
end

function y=normalize(x)
%    y=x;
     s=std(x);
     if s>0.00000001
         y=(x-mean(x))./s;
     else
         y=(x-mean(x));
     end
end
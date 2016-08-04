function [stats,means]=updateStatsAndMean(begends,x,dPointDiffFun,stats,means,motif)
% Updates statistics and mean of a motif (internal function)
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
    x=x(:)';
    n=size(begends,1);
    if n<1
        stats(motif).mean=nan;
        stats(motif).var=nan;
        stats(motif).min=nan;
        stats(motif).max=nan;
        stats(motif).median=nan;
        stats(motif).mode=nan;
        means{motif}=nan;
    elseif n==1
        stats(motif).mean=nan;
        stats(motif).var=nan;
        stats(motif).min=nan;
        stats(motif).max=nan;
        stats(motif).median=nan;
        stats(motif).mode=nan;
        means{motif}=nan;
    end    
    mLen=begends(:,2)-begends(:,1)+1;
    minoLen=min(mLen);  %maxoLen=max(mLen);
    meanMotif=zeros(1,minoLen);    
    for kk=1:n
        meanMotif=meanMotif+(x(begends(kk,1):begends(kk,2)));
    end    
    meanMotif=meanMotif./n;                        
    m=max(begends(:,2)-begends(:,1)+1);
    X=zeros(n,m);       % convert to a cell if the sizes may not be equal
    for ji=1:n
        X(ji,:)=(x(begends(ji,1):begends(ji,2)));
    end
    % @todo here may be we can try to slide the ocurrences to get minimum
    % mean distance overall before calculating the means and variances.
    % his will requuire reading dFun as an input. The same can be done
    % in the other loop calculating newX!!
    d=0; d2=0; l=0;        
    diffs=[];
    for ii=1:n-1
        for ji=ii+1:n                
            diff=dPointDiffFun(X(ii,:)-X(ji,:));
            diffs=[diffs,diff];                
            l=l+length(diff);
            d=d+sum(diff);
            d2=d2+sum(diff.*diff);
        end
    end
    stats(motif).mean=d/l;
    stats(motif).var=(d2-(d^2/l))/(l-1);                    
    stats(motif).min=min(diffs);
    stats(motif).max=max(diffs);
    stats(motif).median=median(diffs);
    stats(motif).mode=mode(diffs);
    means{motif}=meanMotif;
end
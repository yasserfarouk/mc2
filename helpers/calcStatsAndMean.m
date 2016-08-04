function [stats,means]=calcStatsAndMean(begends,x,dPointDiffFun)
% Recalculates the mean and statistics of a given set of candidate motif
% occurrences (internal function)
%
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%
    x=x(:)';
    mLen=begends(:,2)-begends(:,1)+1;
    minoLen=min(mLen);  %maxoLen=max(mLen);
    meanMotif=zeros(1,minoLen);
    for kk=1:size(begends,1)
        meanMotif=meanMotif+x(begends(kk,1)...
            :begends(kk,2));
    end
    n=size(begends,1);
    meanMotif=meanMotif./n;                        
    m=max(begends(:,2)-begends(:,1)+1);
    X=zeros(n,m);       % convert to a cell if the sizes may not be equal
    for ji=1:n
        X(ji,:)=x(begends(ji,1):begends(ji,2));
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
    stats=struct('mean',0,'var',0,'max',0,'min',0,'median',0,'mode',0);
    stats.mean=d/l;
    stats.var=(d2-(d^2/l))/(l-1);                    
    stats.min=min(diffs);
    stats.max=max(diffs);
    stats.median=median(diffs);
    stats.mode=mode(diffs);
    means=meanMotif;
end
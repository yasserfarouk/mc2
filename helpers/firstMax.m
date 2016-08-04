function [bestI,bestV]=firstMax(data,th)
% calculates the first maximum cahnge in every row of the data that is within th
% fraction of the total maximum change of the data
% data: the input data (a single row)
% th: the fraciton of the data maximum that is accepted. default value is
% 0.75
if(nargin<2)
    th=0.75;
end
if(isempty(data))
    bestI=0; bestV=0;
    return
end
if(numel(data)==1)
    bestI=1;
    bestV=data(bestI);
    return
end
diff=data(1,2:size(data,2))-data(1,1:size(data,2)-1);
maxdiff=max(diff);
diff=thin(diff,th*1e-3*maxdiff);
bestI=find(diff>=th*maxdiff,1,'first');
bestV=data(bestI);
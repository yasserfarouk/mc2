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
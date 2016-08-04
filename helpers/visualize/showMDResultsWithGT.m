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

function [figs]=showMDResultsWithGT(x,locsT,locs,useMotifRaising)
% displays results of MD with the ground truth (or two MD algorithms
if ~exist('useMotifRaising','var') || isempty(useMotifRaising)
    useMotifRaising=false;
end

figs=figure;
%maximize(figs);
subplot(2,1,1);
showMDResultsOnData(x,locsT,useMotifRaising);
ylabel('Ground Truth');
subplot(2,1,2);
showMDResultsOnData(x,locs,useMotifRaising);
ylabel('Discovered');
end

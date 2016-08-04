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

function [isover,fraction]=testOverlapping(begend1,begend2)
% checks whether of not two motif stems are overlapping optionally returing
% the fraction of overlap
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
    isover=(begend1(1)>=begend2(1) && begend1(1)<=begend2(2)) || ...
        (begend1(2)>=begend2(1) && begend1(2)<=begend2(2)) || ...
        (begend2(1)>=begend1(1) && begend2(1)<=begend1(2)) || ...
        (begend2(2)>=begend1(1) && begend2(2)<=begend1(2));
    if isover && nargout>1
        fraction=overlapfraction(begend1,begend2);
    elseif nargout>1
        fraction=0;
    end
end
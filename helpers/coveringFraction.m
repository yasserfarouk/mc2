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

function [fraction,reverse,coveredRange]=coveringFraction(covered,covering)
%finds how much of covered is covered by covering (internal function)
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
    if covered(1)>covering(2) || covering(1)>covered(2)
        fraction =0;
        reverse=0;
        coveredRange=[];
        return;
    end
    if (covering(1)<=covered(1) && covering(2)>= covered(2))
        fraction =1;        
        reverse=(covered(2)-covered(1))./(covering(2)-covering(1));        
        coveredRange=covered;
        return;
    end
    coveredRange=[max([covered(1),covering(1)]),min([covered(2),covering(2)])];
    fraction=(coveredRange(2)-coveredRange(1)+1)/(covered(2)-covered(1)+1);
    reverse=(coveredRange(2)-coveredRange(1)+1)/(covering(2)-covering(1)+1);    
end
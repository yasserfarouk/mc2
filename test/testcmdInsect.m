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

d=csvread('./Data/insect15.txt');
%d=x;
c=cpd(d,8,10);
candLocs=findLocsTh(c,0.008);
[locs,means,stats,maxSmallDist,minLargeDist] =sdCMD(c,[400,700],candLocs);
[locsF,meansF,statsF] = mdFinalize(d',[400,700],candLocs,locs,means,stats);

showMDResults(locs,means,stats,2,2);
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

function [foundLocs,trueLocs] = mdq2cpq(locs,locsT)
%Converts the problem of CMD quality measurment into CPD quality
%measurement
%
% INPUTS:
% ======
%   locsT   The true locations of motifs as a cell array of nMotif*1 and
%           each element is a nOccurrences(i)*2 array
%   locs    Found locations in the same form as locsT
%
% OUTPUTS:
% ========
% foundLocs, trueLocs   the boundaries of the motifs in the form needed by
%                       cpquality() function
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

trueLocs=[];
nMotifs=length(locsT);
for i=1:nMotifs
    nO=size(locsT{i},1);
    for j=1:nO
        trueLocs=[trueLocs, locsT{i}(j,:)];
    end
end

found=[];
nMotifs=length(locs);
for i=1:nMotifs
    nO=size(locs{i},1);
    for j=1:nO
        found=[found, locs{i}(j,:)];
    end
end
found=unique(found);
foundLocs=cell(1,1);
foundLocs{1}=found;
end


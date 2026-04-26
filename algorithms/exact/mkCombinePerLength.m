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

function [ locsF,meansF,statsF,dists,maxSmallDist ] = mkCombinePerLength(x, locs,dists)
%Combines the results of MK algorithm runs into sets of motifs
%   Assumes that every pair in the input is of the same length
%   This is true for all MK based algorithms
locsF=cell(0,1);
meansF=cell(0,1);
statsF=cell(0,1);
maxSmallDist=0;
if isempty(locs) || isempty(dists)
    return;
end
[maxSmallDist,locs,dists]=mkRemoveLargeDists(locs,dists);
locsM=cell2mat(locs);
lengths=locsM(:,2)-locsM(:,1)+1;
[lengths,I]=sort(lengths);
locsM=locsM(I,:);
ls=sort(unique(lengths),'descend');
n=numel(ls);

for i=1:n
    thisLength=find(lengths==ls(i));
    dd=dists(thisLength(2:2:end)./2);    
    [lf,mf,sf]=mkCombine(x,locsM(thisLength,:),ls(i),dd);
    locsF=[locsF;lf];
    if nargout>1
        meansF=[meansF,mf];
        if nargout>2
            statsF=[statsF;sf];
        end
    end
end

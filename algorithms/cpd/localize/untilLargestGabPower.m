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

function i=untilLargestGabPower(x,argvarin)
% returns the index after which the largest gab in x appears. x should be
% sorted in descending order. argvarin is not used and is here just as a place holder for the
% function to be the same as all other until* functions
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

n=numel(x);
if(n==0)
    i=0;
    return;
end
if(n==1)
    i=1;
    return;
end
x2=x.*x;
i=1; m=x2(1)-x2(2);
for j=2:n-1
    d=x2(j)-x2(j+1);
    if(d)>m
        m=d;
        i=j;
    end
end
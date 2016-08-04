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

function [x] = randp(p,n)
%RANDP Returns n integers corresponding to samples from p ranging from 1 to
%numel(p)
%   Detailed explanation goes here

P=[cumsum(p)];
P=P./P(end);

r=rand(n,1);
%r=sort(r);
x=zeros(n,1);
x(1)=find(P>=r(1),1);
for i=2:n
    %x(i)=x(i-1)+find(P(x(i-1):end)>=r(i),1);
    x(i)=find(P>=r(i),1);
end
%x=x(randperm(n));

end


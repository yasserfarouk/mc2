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

figure; 
aa=[
    -0.7,0.5,0.9;
    -0.2,0.5,0.9;
    -0.7,1.5,0.09;
    -0.0,0.5,0.04;
    ];
for i=1:4; 
    subplot(2,2,i); 
    a=aa(i,:)'; 
    x=generateAR(a,100,1); plot(x); 
    title(sprintf('%3.2f, ',a));  
end;
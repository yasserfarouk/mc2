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

A=[2.7607 -3.8106 2.6535 -0.9238];

rng default;
x=generateAR(A,1024,1,[],0.04,0); 

[a,sigma2,x_0_m,a0,mu_epsilon]=learnARLS(x,4); a=a';
%[ay,sigma2y,x_0_my,a0y,mu_epsilony]=learnARYuleWalker(x,4);
[ay,sigma2y]=aryule(x,4); ay=-ay(2:end);
rng default;
y=generateAR(a,1024,1,[],sigma2,0);

rng default;
z=generateAR(ay,1024,1,[],sigma2y,0);

figure;
subplot(2,1,1);
plot([x,y,z]);
legend({'original','least squares','Yule-Walker'});
title('original signal and learned predictions');
subplot(2,1,2);
plot([y-x,z-x]);
legend({'least squares','Yule-Walker'});
title('Learning error');
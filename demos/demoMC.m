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
subplot(2,2,1); 
plot(generateMarkovChain([0;0],eye(2),eye(2),100),'LineWidth',2);
title('$\Sigma=I$','interpreter','latex', 'FontSize',20);
subplot(2,2,2); 
plot(generateMarkovChain([0;0],eye(2),[0.2,0.0;0.0,0.6],100),'LineWidth',2);
title('$\Sigma  = \left[ {\begin{array}{*{20}{c}}{0.2}&0 \\ 0&{0.6}\end{array}} \right]$','interpreter','latex','FontSize',20);
subplot(2,2,3); 
plot(generateMarkovChain([0;0],eye(2),[0.2,0.1;0.1,0.6],100),'LineWidth',2);
title('$\Sigma  = \left[ {\begin{array}{*{20}{c}}{0.2}&0.1 \\ 0.1&{0.6}\end{array}} \right]$','interpreter','latex','FontSize',20);

subplot(2,2,4); 
plot(generateMarkovChain([0;0],eye(2),[0.2,0.3;0.3,0.6],100),'LineWidth',2);
title('$\Sigma  = \left[ {\begin{array}{*{20}{c}}{0.2}&0.3 \\ 0.3 &{0.6}\end{array}} \right]$','interpreter','latex','FontSize',20);

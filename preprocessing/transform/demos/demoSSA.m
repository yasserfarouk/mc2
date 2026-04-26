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

groundTruth=sin([1:500]*pi/50)+2*sin([1:500]*pi/5)+[-249:250]./50;
noise=randn(1,500);
x=groundTruth+noise;

[h,zeta]=ssaDecompose(x,50);

figure;
subplot(4,2,1); plot(groundTruth,'LineWidth',2); title('Ground Truth', 'FontSize',20);
subplot(4,2,2); plot(x,'LineWidth',2); title('Input after adding Gaussian noise', 'FontSize',20);
for i=1:6
    subplot(4,2,i+2);
    plot(hankel2ts(h{i}),'LineWidth',2);
    title(sprintf('expansion time-series #%d (weight=%2.3f)',i,zeta(i)),'FontSize',20);
end

[x2,zeta]=ssaReconstruct(h,zeta,0.4);
figure;
subplot(2,2,1); plot(groundTruth,'LineWidth',2); title('Ground Truth', 'FontSize',20);
subplot(2,2,2); plot(x,'LineWidth',2); title('Input after adding Gaussian noise', 'FontSize',20);
subplot(2,2,3); plot(x2{1},'LineWidth',2); title('Approximation using SSA ','FontSize',20);
subplot(2,2,4); plot(x2{2},'LineWidth',2);
hold on;
plot(noise,'r','LineWidth',2);
title('Estimated (blue) and real (red) noise)', 'FontSize',20);
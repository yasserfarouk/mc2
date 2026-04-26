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

groundTruth=generateMarkovChain(zeros(1,2), eye(2),[0.2,0.1;0.1,0.6],500)';
noise=randn(2,500);
x=groundTruth+noise;

% using SSAM
[h,zeta]=ssaDecompose(x,50);

figure;
subplot(4,2,1); plot(groundTruth','LineWidth',2); title('Ground Truth', 'FontSize',20);
subplot(4,2,2); plot(x','LineWidth',2); title('Input after adding Gaussian noise', 'FontSize',20);
for i=1:6
    subplot(4,2,i+2);
    plot(hankel2ts(h{i},2)','LineWidth',2);
    title(sprintf('expansion time-series #%d (weight=%2.3f)',i,zeta(i)),'FontSize',20);
end

[x2,zeta]=ssaReconstruct(h,zeta,0.75,2);
fgg=figure;
subplot(2,2,1); plot(groundTruth','LineWidth',2); title('Ground Truth', 'FontSize',20);
subplot(2,2,2); plot(x','LineWidth',2); title('Input after adding Gaussian noise', 'FontSize',20);
subplot(2,2,3); plot(x2{1}','LineWidth',2); title('SSAM','FontSize',20);
%subplot(2,2,4); plot(x2{2}','LineWidth',2);
hold on;
plot(noise,'r','LineWidth',2);
title('Estimated (blue) and real (red) noise)', 'FontSize',20);

% Repeated use of SSA
X=x;
GT=groundTruth;
fg1=figure;
x4=x;
for m=1:2
    x=X(m,:);
    groundTruth=GT(m,:);
    [h,zeta]=ssaDecompose(x,50);
    figure(fg1);
    subplot(4,2,1); hold on; plot(groundTruth,'LineWidth',2); title('Ground Truth', 'FontSize',20);
    subplot(4,2,2); hold on; plot(x,'LineWidth',2); title('Input after adding Gaussian noise', 'FontSize',20);
    for i=1:6
        subplot(4,2,i+2);
         hold on;
        plot(hankel2ts(h{i}),'LineWidth',2);
        title(sprintf('expansion time-series #%d (weight=%2.3f)',i,zeta(i)),'FontSize',20);
    end

    [x3,zeta]=ssaReconstruct(h,zeta,0.75);
    figure(fgg);
   
    subplot(2,2,4);  hold on; plot(x3{1},'LineWidth',2); title('Approximation using SSA ','FontSize',20);
    %subplot(2,2,4);  hold on; plot(x3{2},'LineWidth',2);
    hold on;
    plot(noise,'r','LineWidth',2);
    title('Repeated SSA', 'FontSize',20);
    
    x4(m,:)=x3{1};
end
fprintf('Error in input: %5.3f\n', norm(GT-X));
fprintf('Error in SSAM approximation: %5.3f\n',norm(GT-x2{1}));
fprintf('Error in repeated SSA: %5.3f\n',norm(GT-x4));
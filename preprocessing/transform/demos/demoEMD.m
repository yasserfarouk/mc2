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

imf=emd(x);
n=size(imf,1);
nRow=ceil(n/2)+1;
nCol=2;
figure;
subplot(nRow,nCol,1); plot(groundTruth,'LineWidth',2); title('Ground Truth', 'FontSize',20);
subplot(nRow,nCol,2); plot(x,'LineWidth',2); title('Input after adding Gaussian noise', 'FontSize',20);
for i=1:n
    subplot(nRow,nCol,i+2);
    plot(imf(i,:),'LineWidth',2);
    title(sprintf('IMF #%d',i),'FontSize',20);
end




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

N=200;
L=100;
T=3000;
range=[-5,5];    
scale2range=1;
noiseLevel=0.00;
gNoiseSigmaInitial=1.0;       % initial noise
pOutlier=0.0;
minChangeDistance=150;
pchange=0.01;


[xbase,locsTbase,tbase,xTruebase]=produceSingle(T,'pchange',pchange,...
    'minChangeDistance',minChangeDistance,'noiseLevel',noiseLevel,...
    'poutlier',pOutlier,'noiseSigma',gNoiseSigmaInitial,'range',range...
    ,'scale2range',scale2range);
xbase=[xbase,randn(1,numel(xbase))];
tbase=[tbase,zeros(1,numel(tbase))];
y=sst(xbase,10,50,3,0,0);
z=rsst(xbase,10,50);

figure;
subplot(3,1,1);
plot(xbase);
hold on;
plot(tbase);
ylabel('$x_t$');

subplot(3,1,2);
plot(y./max(y));
hold on;
plot(tbase);
ylabel('SST');

subplot(3,1,3);
plot(z./max(z));
hold on;
plot(tbase);
ylabel('RSST');
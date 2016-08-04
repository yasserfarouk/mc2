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

l=0.2;
sigmaN=0.1;
predictWithNoise=true;
const=1;
backcolor=[0.9,0.9,0.9];
%fMean=@(x) (meanSinusoidal(x,[const;4*ones(size(x,1),1);4*zeros(size(x,1),1)]));
fMean=@(x) (meanLinear(x));
X=-2:0.01:1;
n=numel(X);
M=10;
f=zeros(M,n);
for m=1:M
    [f(m,:),K]=generateGP(X,sigmaN,@covSE,l,fMean);
end
figure;
subplot(2,2,1);
hold on;
show1DGaussianFun(fMean(X),K,X,[0,0,0],backcolor);
plot(X,f);
xlabel('(a) Before seeing any inputs','FontSize',16);

%Xin = [-1.5 -1 -0.75 -0.4 -0.3 0];
%f = [-1.6 -1.3 -0.5 0 0.3 0.6];

Xin=[-1.5,-0.3,0.38];
f=[-1.6,0.3,-4];

gp=createGP([],[],sigmaN,@covSE,l,fMean);
chr={'b','c','d'};
for i=1:3
    gp=add2GP(gp,Xin(i),f(i));
    Xstar=X;
    fstar=zeros(M,n);
    [mu,K]=generateFromGP(gp,Xstar,predictWithNoise);
    for m=1:M
        [fstar(m,:)]=grand(mu,K)';
    end
    subplot(2,2,i+1);
    hold on;
    show1DGaussianFun(mu,K,Xstar,[0,0,0],backcolor);
    plot(X,fstar);
    plot(Xin(1:i),f(1:i),'+r','LineWidth',4);
    xlabel(sprintf('(%s) Input %d: log p(y|X)=%4.3f',chr{i},i,gp.logpyX),'FontSize',16);
end
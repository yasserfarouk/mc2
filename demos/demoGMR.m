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
K=4;
N=4;
prior=[0.2470    0.2366    0.2512    0.2652];
mu=[1.7921    1.8765    1.8903    1.9371;
    2.1489    1.9193    2.0379    1.9056;
    1.6503    1.7524    1.8556    1.8100;
    2.2006    1.9675    2.0135    1.7831;
    6.3562   41.9089   30.4847   18.3157];
sigma=zeros(5,5,4);
sigma(:,:,1)=[0.0096   -0.0120    0.0066   -0.0111    0.1046;
   -0.0120    0.0188   -0.0082    0.0166   -0.2195;
    0.0066   -0.0082    0.0063   -0.0076   -0.0190;
   -0.0111    0.0166   -0.0076    0.0156   -0.1958;
    0.1046   -0.2195   -0.0190   -0.1958   11.7227
];
sigma(:,:,2) =[
    0.0157   -0.0298    0.0108   -0.0270    0.3160;
   -0.0298    0.0824   -0.0265    0.0749   -0.7203;
    0.0108   -0.0265    0.0127   -0.0297    0.2893;
   -0.0270    0.0749   -0.0297    0.0817   -0.8342;
    0.3160   -0.7203    0.2893   -0.8342   10.5598];


sigma(:,:,3) =[
    0.0219   -0.0298    0.0241   -0.0342   -0.3789;
   -0.0298    0.0486   -0.0371    0.0532    0.6325;
    0.0241   -0.0371    0.0314   -0.0452   -0.5482;
   -0.0342    0.0532   -0.0452    0.0686    0.8202;
   -0.3789    0.6325   -0.5482    0.8202   12.1554
];

sigma(:,:,4) =[
    0.0075   -0.0173    0.0030   -0.0066   -0.2085;
   -0.0173    0.0507   -0.0094    0.0225    0.5717;
    0.0030   -0.0094    0.0047   -0.0058   -0.0694;
   -0.0066    0.0225   -0.0058    0.0147    0.1848;
   -0.2085    0.5717   -0.0694    0.1848   13.7319
];

letters={'(a)','(b)','(c)','(d)'};
 
plot(generateGMR(prior,mu,sigma,47),'LineWidth',2);

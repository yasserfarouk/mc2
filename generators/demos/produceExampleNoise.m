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

randPP=0.8;
M=[2,3,6;3,5,5;4,6,11;5,1,4;3,7,10];
[X,Y,P,t,d]=produce(10,30,M); X=X-0.5*randPP+randPP*rand(size(X)); Y=Y-0.5*randPP+randPP*rand(size(Y));  subplot(2,1,1); plot(X); subplot(2,1,2); plot(Y); 
savefig('./Data/exampleN','png', '-rgb', '-c0.1', '-r250');
%copyfile('exampleN.png','C:\Users\yasser\Research\Papers\Mypapers\_current\IntelligentSystemsKnowldgeBook\exampleN.png');
E=zeros(size(M,1),3);
E(:,1)=1:size(M,1);
E(:,2)=M(:,1);
E(:,3)=M(:,2);
grPlot([],E,'d','%d','%d');
savefig('./Data/exampleGraphN','png', '-rgb', '-c0.1', '-r250');
%copyfile('exampleGraphN.png','C:\Users\yasser\Research\Papers\Mypapers\_current\IntelligentSystemsKnowldgeBook\exampleGraphN.png');
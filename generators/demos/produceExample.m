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

M=[2,3,6;3,5,5;4,6,11;5,1,4;3,7,10];
[X,Y,P,t,d]=produce(100,30,M); subplot(2,1,1); plot(X); subplot(2,1,2); plot(Y); 
savefig('example','png', '-rgb', '-c0.1', '-r250');

E=zeros(size(M,1),3);
E(:,1)=1:size(M,1);
E(:,2)=M(:,1);
E(:,3)=M(:,2);
grPlot([],E,'d','%d','%d');
savefig('./Figs/exampleGraph','png', '-rgb', '-c0.1', '-r250');
%copyfile('exampleGraph.png','C:\Users\yasser\Research\Papers\Mypapers\_current\IntelligentSystemsKnowldgeBook\exampleGraph.png');
%copyfile('example.png','C:\Users\yasser\Research\Papers\Mypapers\_current\IntelligentSystemsKnowldgeBook\example.png');
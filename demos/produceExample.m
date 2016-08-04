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
function [x] = randp(p,n)
%RANDP Returns n integers corresponding to samples from p ranging from 1 to
%numel(p)
%   Detailed explanation goes here

P=[cumsum(p)];
P=P./P(end);

r=rand(n,1);
%r=sort(r);
x=zeros(n,1);
x(1)=find(P>=r(1),1);
for i=2:n
    %x(i)=x(i-1)+find(P(x(i-1):end)>=r(i),1);
    x(i)=find(P>=r(i),1);
end
%x=x(randperm(n));

end


function x=generateMA(a,T,n)
% generates data from a moving average process of order m [MA(m)]
%
% a     The parameters of the MA(m) process and m=numel(a) in the order
%       [a_0, a_1, ... a_m]
% T     number of poitns in the time series
% n     number of time series (default 1)
%
%
% output:
% =======
% x     T*n matrix representing n different time-series each of length n
%       generated from a standard Gaussian white noise

if ~exist('n','var') || isempty(n)
    n=1;
end
m=numel(a)-1;
theta=randn(T+m,n);
x=zeros(T,n);
a=a(:);
for i=1:T
    x(i,:)=sum(theta(i+m:-1:i,:).*repmat(a,1,n));
end

end
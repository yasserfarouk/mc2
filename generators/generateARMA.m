function [x,noise]=generateARMA(a,b,T,n,sigma2,x_0_m,a0,ignoreFirstM)
% generates data from an ARMA(p,n) process consisting of AR(p) and MA(n)
% processes
%
% a     The parameters of the AR(p) process and p=numel(a) in the order
%       [a_1, ... a_m]
% b     The parameters of the MA(n) process and p=numel(a) in the order
%       [b_0, ... b_m]
% T     number of poitns in the time series
% n     number of time series (default 1)
%
% x_0_m p*n matrix giving the initial values of the first p
%       elements of  each time-series. If not given they are selected from a standard
%       normal distribution. If a vector is given it is then used to
%       initialize all of the n time-series
% sigma2 The variance of white noise to be added to the time-series.
%        Default is 1. If zero no white noise is added
%
% a_0    A constant value to be added to the time-series. 
% 
% ignoreFirstM  if nonzero then the first p samples are created then ignored
%
% output:
% =======
% x     T*n matrix representing n different time-series each of length n
%       generated from a standard Gaussian white noise
%
%
% Notes:
% ======
% The first p values are assumed to be zeros.
p=numel(a);
q=numel(b);
if ~exist('n','var') || isempty(n)
    n=1;
end
if ~exist('sigma2','var') || isempty(sigma2)
    sigma2=1.0;
end
sigma=sqrt(sigma2);

if ~exist('a0','var') || isempty(a0)
    a0=0.0;
end
if numel(a0)==1
    a0=repmat(a0,1,n);
end
noise=zeros(p,n);
rmin=min([p,q]);
rmax=max([p,q]);
if ~exist('x_0_m','var') || isempty(x_0_m)
    x_0_m=sigma.*randn(rmax,n);
    noise=x_0_m;
    for i=2:rmin
        x_0_m(i,:)=a0+x_0_m(i,:)+sum(x_0_m(i-1:-1:1,:).*...
            repmat(a(1:i-1)',1,n)) ...
            +sum(noise(i:-1:1,:).*...
            repmat(b(1:i)',1,n));
    end
    for i=rmin+1:q
        x_0_m(i,:)=a0+x_0_m(i,:)+ ...
            +sum(noise(i:-1:1,:).*...
            repmat(b(1:i)',1,n));
    end
    for i=rmin+1:p
        x_0_m(i,:)=a0+x_0_m(i,:)+sum(x_0_m(i-1:-1:1,:).*...
            repmat(a(1:i-1)',1,n));
    end
end
if ~exist('ignoreFirstM','var') || isempty(ignoreFirstM)
   ignoreFirstM=false;
end
if isvector(x_0_m)    
    x_0_m=repmat(x_0_m(:),1,n);
end
if sigma>1e-6
    x=[x_0_m;sigma.*randn(n,T)'];
else
    x=[x_0_m;zeros(T,n)];
end
noise=[noise;x(rmax+1:end,:)];
a=a(:); b=b(:);
for i=rmax+1:T+rmax
    x(i,:)=a0+sum(x(i-1:-1:i-p,:).*repmat(a,1,n))+...
        sum(noise(i:-1:i-q+1,:).*repmat(b,1,n));
end

if ignoreFirstM
    x=x(rmax+1:end,:);
    noise=noise(rmax+1:end,:);
else
    x=x(1:T,:);
    noise=noise(1:T,:);
end
end
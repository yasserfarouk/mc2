function x=generateARIMA(a,b,d,T,n,sigma2,x_0_m,a0,ignoreFirstM)
% generates data from an ARMA(m,n) process consisting of AR(m) and MA(n)
% processes
%
% a     The parameters of the AR(m) process and m=numel(a) in the order
%       [a_1, ... a_m]
% b     The parameters of the MA(n) process and m=numel(a) in the order
%       [b_0, ... b_m]
% d     The number of differentiations to get the underlying ARMA model
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
% output:
% =======
% x     T*n matrix representing n different time-series each of length n
%       generated from a standard Gaussian white noise
%
%
% Notes:
% ======
% The first m values are assumed to be zeros.
m=numel(a);
if ~exist('n','var') || isempty(n)
    n=1;
end
if ~exist('x_0_m','var') || isempty(x_0_m)
    x_0_m=[];
end
if ~exist('sigma2','var') || isempty(sigma2)
    sigma2=[];
end

if ~exist('a0','var') || isempty(a0)
    a0=[];
end
if ~exist('ignoreFirstM','var') || isempty(ignoreFirstM)
   ignoreFirstM=[];
end
method=2;

switch(method)
    case 1
        x=generateARMA(a,b,T,n,sigma2,x_0_m,a0,ignoreFirstM);
        for i=1:d
            x=cumsum(x);
        end
    case 2
        for i=1:d
            a=[a,0]-[1,a];
        end
        x=generateARMA(a,b,T,n,sigma2,x_0_m,a0,ignoreFirstM);
end
end
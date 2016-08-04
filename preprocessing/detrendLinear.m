function [y,a,b]=detrendLinear(x,t)
% removes a linear trend from x
%
% x     A vector representing a time-series
% t     The independent variable. Must be a vector of the same length as x

% y     The time-sereis after detrending
% a,b   parameters of the line used for detrending (x=at+b)

assert(isvector(x),'x must be a vector');

T=numel(x);
if ~exist('t','var') ||  isempty(t)
    t=1:T;
end
B=x(:);
t=t(:);
A=[t,ones(T,1)];
line=A\B;
a=line(1);
b=line(2);
z=a.*t+b;
y=B-z;
if size(x,1)==1
    y=y';
end
end





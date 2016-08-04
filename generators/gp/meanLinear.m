function [ f ] = meanLinear( X,hyp )
%MEANLINEAR Summary of this function goes here
%   Detailed explanation goes here
if isempty(X)
    f=[];
    return;
end
D=size(X,1);
n=size(X,2);
if nargin<2|| isempty(hyp)
    hyp=zeros(D+1,1);
end
w=hyp(2:end);
w0=hyp(1);

f=w0.*ones(1,n)+(w*X')';

end


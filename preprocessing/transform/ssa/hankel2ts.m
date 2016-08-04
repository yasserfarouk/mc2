function x=hankel2ts(H,M)
if nargin<2
    M=1;
end
[~,x]=hankelize(H,M);
end
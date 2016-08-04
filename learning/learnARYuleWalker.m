function [a,sigma2,x_0_m,a0,mu_epsilon]=learnARYuleWalker(x,p,removeMean)
% learns the parameter vector a and covariance of Gaussian noise generating
% an input time-series
if ~exist('removeMean','var')
    removeMean=false;
end
[T,n]=size(x);
x_0_m=x(1:p,:);
if removeMean
    a0=mean(x,1);
    x=x-repmat(a0,T,1);
else
    a0=zeros(1,n);
end
nr=2*p;  % only half of these need to be really calculated because of symmetry
rxx=zeros(nr,1);
for i=1:nr
    rxx(i)=pcf(x,i-p);
end

X=zeros(p,p);
B=rxx(end-p+1:end);
for i=1:p
    X(i,:)=rxx(i:i+p-1);
end

a=X\B;

eps=zeros(n*(T-p+1),1);
for i=1:n
    for j=p+1:T
        eps(j,i)=x(j,i)-sum(a.*x(j-p:j-1,i));
    end
end
a=a(end:-1:1);
sigma2=var(eps(:));
if nargout>3
    mu_epsilon=mean(eps(:));
end


a=a(:)';
end

function rxx_k=pcf(x,k)
    if k>-1
        xx=x(k+1:end,:).*x(1:end-k,:);
    else
        xx=x(1-k:end,:).*x(1:end+k,:);
    end
    xx=xx(:);
    rxx_k=mean(xx);
end
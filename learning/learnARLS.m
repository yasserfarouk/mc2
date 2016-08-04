function [a,sigma2,x_0_m,a0,mu_epsilon,eps]=learnARLS(x,p,removeMean)
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
X=zeros(n*(T-p+1),p);
B=zeros(n*(T-p+1),1);
nxt=1;
for i=1:n
    for j=p+1:T
        X(nxt,:)=x(j-p:j-1,i);
        B(nxt)=x(j,i);
        nxt=nxt+1;
    end
end

a=X\B;
eps=[x_0_m;zeros(n*(T-p+1),1)];
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
function [a,b,sigma2,x_0_m,a0,mu_epsilon,eps]=learnARMALS(x,p,q,removeMean)
% learns the parameter vector a and covariance of Gaussian noise generating
% an input time-series
[T,n]=size(x);
assert(q>1,'q must be at least two otherwise you can use an AR modelS');
if ~exist('removeMean','var')
    removeMean=false;
end

if removeMean
    a0=mean(x,1);
    x=x-repmat(a0,T,1);
end

[a,sigma2,x_0_m,a0,mu_epsilon,eps]=learnARLS(x,p+q,false);

r=p+q;
rmax=max([p,q]);
X=zeros(n*(T-p+1),r);
B=zeros(n*(T-p+1),1);
nRep=1;

for rep=1:nRep

    nxt=1;
    for i=1:n
        for j=rmax+1:T
            X(nxt,:)=[x(j-1:-1:j-p,i);eps(j-1:-1:j-q,i)];
            B(nxt)=x(j,i)-eps(j,i);
            nxt=nxt+1;
        end
    end

    a=X\B;
    b=a(p+1:end);
    a=a(1:p);
    %eps=zeros(n*(T-p+1),1);
    for i=1:n
        for j=rmax+1:T
            eps(j,i)=x(j,i)-sum(a.*x(j-1:-1:j-p,i))-sum(b(2:end).*eps(j-1:-1:j-q+1,i));
        end
    end
    
end;

sigma2=var(eps);
if nargout>3
    mu_epsilon=mean(eps);
end
a=a(:)';
b=b(:)';

end
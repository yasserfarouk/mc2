function [p0,A,mu,sigma,log_pX,err,iter]=learnHMM(x,K,MaxIter,MinErr,p0,A,mu,sigma)
% learns the parameters of a HMM from an input time-series using Baum Welch
% algorithm

[T,n]=size(x);
if ~exist('MaxIter','var') || isempty(MaxIter)
    MaxIter=1000;
end
if ~exist('MinErr','var') || isempty(MinErr)
    MinErr=1e-4;
end
[idx,C]=kmeans(x,K,'EmptyAction','singleton');
if ~exist('p0','var') || isempty(p0)
    p0=ones(K,1)./T;
    for i=1:K
        p0(i)=p0(i).*numel(find(idx==i));
    end
end

if ~exist('mu','var') || isempty(mu)
    %mu=x(randi(T,K,1),:)';
    mu=C';
end

if ~exist('sigma','var') || isempty(sigma)
    sigma=zeros(n,n,K);
    for k=1:K
        data=x(idx==k,:);
        nd=size(data,1);
        for i=1:n
            di=data(:,i)-mu(i,k);
            for j=1:n
                dj=data(:,j)-mu(j,k);
                sigma(i,j,k)= sum((di.*dj)) /(nd);
            end
        end
    end
end
if ~exist('A','var') || isempty(A)
    A=zeros(K,K);
    for i=1:T-1
        A(idx(i),idx(i+1))=A(idx(i),idx(i+1))+1;
    end
    for i=1:K
        A(i,:)=A(i,:)./sum(A(i,:));
    end
end
smallNum=1e-6;
for iter=1:MaxIter
    alpha=zeros(K,T);
    pObservations=zeros(K,T);
    for k=1:K
        pObservations(k,:)=gaussPDF(x,mu(:,k),sigma(:,:,k));
    end
    c=zeros(T,1);
    % forward recursion
    alpha(:,1)=p0.*pObservations(:,1);
    c(1)=sum(alpha(:,1));
    alpha(:,1)=alpha(:,1)./c(1);
    for t=1:T-1
        for k=1:K
            alpha(k,t+1)=pObservations(k,t+1).*(A(:,k)'*alpha(:,t))+smallNum;
        end
        c(t+1)=sum(alpha(:,t+1));
        alpha(:,t+1)=alpha(:,t+1)./c(t+1);
    end
    log_pX=-(sum(log(c)));
    % backward recursion
    beta=zeros(K,T);
    beta(:,T)=ones(K,1);

    for t=T-1:-1:1
        for k=1:K
            beta(k,t)=(A(k,:)*(pObservations(:,t+1).*beta(:,t+1)))+smallNum;
        end
        beta(:,t)=beta(:,t)./sum(beta(:,t));
    end

    gamma=alpha.*beta;
    for t=1:T
        gamma(:,t)=gamma(:,t)./sum(gamma(:,t));
    end

    zeta=zeros(K,K,T-1);
    for t=1:T-1
        for i=1:K
            for j=1:K
                zeta(i,j,t)=alpha(i,t).*A(i,j).*pObservations(j,t+1).*beta(j,t+1);
            end
        end
        zeta(:,:,t)=zeta(:,:,t)./sum(sum(zeta(:,:,t)));
    end
    p0old=p0;
    muold=mu;
    sigmaold=sigma;
    Aold=A;

    p0=gamma(:,1);
    for i=1:K
        %s=sum(gamma(i,1:end-1));
        for j=1:K
            A(i,j)=sum(zeta(i,j,:));
        end
    end
    for i=1:K
        A(i,:)=A(i,:)./sum(A(i,:));
    end
    for k=1:K
        sg=sum(gamma(k,:));
        for i=1:n
            mu(i,k)=(gamma(k,:)*x(:,i))./sg;
        end
    end

    for k=1:K
        sg=sum(gamma(k,:));
        diff=(x-repmat(mu(:,k)',T,1))';
        sigma(:,:,k)=zeros(n,n);
        for t=1:T
            ss=diff(:,t)*diff(:,t)';
            sigma(:,:,k)= sigma(:,:,k)+gamma(k,t).*(ss);
        end
        sigma(:,:,k)=sigma(:,:,k)./sg;
    end

    err=(norm(abs(Aold-A))+sum(abs(p0old-p0))+sum(sum(abs(muold-mu)))+sum(sum(sum(abs(sigmaold-sigma)))))./(K+n*K+n*n*K+K*K);
    if err<MinErr
        break;
    end
end
end

function p=gaussPDF(x,mu,sigma)
    [T,n] = size(x);
    x = x - repmat(mu',T,1);
    p = sum((x*inv(sigma)).*x, 2);
    p = exp(-0.5*p) / sqrt((2*pi)^n * (abs(det(sigma))+realmin));
end
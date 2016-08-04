function [p0,A,mu,sigma,pX,err,iter]=learnHMMMulti(xs,K,combinationMethod,MaxIter,MinErr,p0,A,mu,sigma)
% learns the parameters of a HMM from an input time-series using Baum Welch
% algorithm
assert(iscell(xs),'xs must be a cell vector with each time-series in a cell and all time-series must have the same number of dimensions (e.g. size(x{i},2)=size(x{j},2) for any i,j)');
nSeries=numel(xs);
xAll=[];
for s=1:nSeries
    xAll=[xAll;xs{s}];
end
if ~exist('combinationMethod','var') || isempty('combinationMethod')
    combinationMethod='average';
    combinationMethod='weighted-by-1/p';
    combinationMethod='weighted-by-p';
end
if ~exist('MaxIter','var') || isempty(MaxIter)
    MaxIter=1000;
end
if ~exist('MinErr','var') || isempty(MinErr)
    MinErr=1e-4;
end
[T,n]=size(xAll);
[idx,C]=kmeans(xAll,K,'EmptyAction','singleton');
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
        data=xAll(idx==k,:);
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
p0s=zeros(K,nSeries);
As=zeros(K,K,nSeries);
mus=zeros(n,K,nSeries);
sigmas=zeros(n,n,K,nSeries);
log_pXs=zeros(nSeries,1);
for iter=1:MaxIter
    p0old=p0;
    muold=mu;
    sigmaold=sigma;
    Aold=A;
    gamma=cell(nSeries,1);
    zeta=cell(nSeries,1);
    for s=1:nSeries
        x=xs{s};
        [T,n]=size(x);
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
        log_pXs(s)=-(sum(-log(c)));
        % backward recursion
        beta=zeros(K,T);
        beta(:,T)=ones(K,1);

        for t=T-1:-1:1
            for k=1:K
                beta(k,t)=(A(k,:)*(pObservations(:,t+1).*beta(:,t+1)))+smallNum;
            end
            beta(:,t)=beta(:,t)./sum(beta(:,t));
        end

        gamma{s}=alpha.*beta;
        for t=1:T
            gamma{s}(:,t)=gamma{s}(:,t)./sum(gamma{s}(:,t));
        end

        zeta{s}=zeros(K,K,T-1);
        for t=1:T-1
            for i=1:K
                for j=1:K
                    zeta{s}(i,j,t)=alpha(i,t).*A(i,j).*pObservations(j,t+1).*beta(j,t+1);
                end
            end
            zeta{s}(:,:,t)=zeta{s}(:,:,t)./sum(sum(zeta{s}(:,:,t)));
        end
        
    end
    switch(combinationMethod)
        case 'average'
            pX=ones(nSeries,1)./nSeries;
        case 'weighted-by-1/p'
            pXsum=elogsum(log_pXs);
            pX=log_pXs-pXsum;
            pX=1/(exp(pX)+realmin);
            pX=pX./sum(pX);
        case 'weighted-by-p'
            pXsum=elogsum(log_pXs);
            pX=log_pXs-pXsum;
            pX=(exp(pX)+realmin);
            pX=pX./sum(pX);
        otherwise
    end
    
    %pX=zeros(size(pX)); pX(1)=1;
    
    p0=zeros(size(p0));
    for s=1:nSeries
        p0=p0+gamma{s}(:,1).*pX(s);
    end
    p0=p0./sum(p0);
    
    A=zeros(size(A));
    for s=1:nSeries
        for i=1:K
            %s=sum(gamma(i,1:end-1));
            for j=1:K
                A(i,j)=A(i,j)+sum(zeta{s}(i,j,:)).*pX(s);
            end
        end
    end
    for i=1:K
        A(i,:)=A(i,:)./sum(A(i,:));
    end
    
    mus=cell(nSeries,1);
    for s=1:nSeries
        mus{s}=mu;
        for k=1:K
            sg=sum(gamma{s}(k,:));
            for i=1:n
                mus{s}(i,k)=((gamma{s}(k,:)*(xs{s}(:,i)))./sg);
            end
        end
    end
    mu=zeros(size(mu));
    for s=1:nSeries
        mu=mu+mus{s}.*pX(s);
    end
    
    sigmas=cell(nSeries,1);
    for s=1:nSeries
        sigmas{s}=zeros(size(sigma));
        T=size(xs{s},1);
        for k=1:K
            sg=sum(gamma{s}(k,:));
            diff=(xs{s}-repmat(mu(:,k)',T,1))';
            %sigma(:,:,k)=zeros(n,n);
            for t=1:T
                ss=diff(:,t)*diff(:,t)';
                sigmas{s}(:,:,k)= sigmas{s}(:,:,k)+(gamma{s}(k,t).*(ss));
            end
            sigmas{s}(:,:,k)=sigmas{s}(:,:,k)./sg;
        end
    end
    sigma=zeros(size(sigma));
    for s=1:nSeries
        sigma=sigma+sigmas{s}.*pX(s);
    end
    for k=1:K
        sigma(:,:,k)=nearestSPD(sigma(:,:,k));
    end
    err=(sum(abs(p0old-p0))+sum(sum(abs(muold-mu)))+sum(sum(sum(abs(sigmaold-sigma)))))./(K+n*K+n*n*K);
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

function y=eexp(x)
    if isnan(x) 
        y=0;
    else
        y=exp(x);
    end
end
function y=elog(x)
    if x==0
        y=nan;
    elseif x>0 
        y=log(x);
    else
        error('elog of a negative number');
    end
end
function lgz=elogsum(lgx,lgy)
    if nargin==1 && isvector(lgx)
        lgz=lgx(1);
        for i=2:numel(lgx)
            lgz=elogsum(lgz,lgx(i));
        end
    else
        if isnan(lgx)
            lgz=lgy;
        elseif isnan(lgy)
            lgz=lgx;
        else
            if lgx>lgy
                lgz=lgx+elog(1+exp(lgy-lgx));
            else
                lgz=lgy+elog(1+exp(lgx-lgy));
            end
        end
    end
end
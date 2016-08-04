function [x,s]=generateHMM(A,p0,mu,sigma,T,a)
% generates data from a HMM using p0 as initial probability and A
% as the transition matrix and Gaussian observation probabilities
%
% A         a square Ns*Ns matrix giving the transition probability from row state to column state    
% p0        an Ns elements vector initial probability distribution
% mu        An N*Ns matrix where mu(:,i) gives the mean of the observation Gaussian for state i
% sigma     An N*N*Ns matrix giving the covariances of observation
%           Gaussians
% T         The lenght of the output time-series
% a     [optional] The parameters of the MA(m) process and m=numel(a) in the order
%       [a_0, a_1, ... a_m]
if ~exist('a','var') || isempty(a)
    a=[];
end
Ns=size(A,1);
assert(size(A,2)==Ns);
assert(all(abs(sum(A,2)-ones(Ns,1))<1e-6));
if ~exist('p0','var') || isempty(p0)
    p0=ones(Ns,1)/Ns;
end
n=size(mu,1);
s=zeros(T,1);
if isempty(a)
    ma=zeros(T,n);
else
    ma=generateMA(a,T,n);
end
if n==1
    x=zeros(T,1);    
    r=randn(T,1);
    s(1)=sample(p0);
    x(1)=mu(s(1))+sqrt(sigma(s(1)))*r(1)+ma(1);
    for t=2:T
        s(t)=sample(A(s(t-1),:));
        x(t)=mu(s(t))+sqrt(sigma(s(t)))*r(t)+ma(t);
    end
else
    x=zeros(T,n);
    s(1)=sample(p0);
    x(1,:)=mvnrnd(mu(:,s(1)),(sigma(:,:,s(1))))+ma(1,:);
    for t=2:T
        s(t)=sample(A(s(t-1),:));
        x(t,:)=mvnrnd(mu(:,s(t)),sigma(:,:,s(t)))+ma(t,:);
    end
end

end

function x=sample(p)
    c=cumsum(p);
    s=c(end)*rand(1,1);
    x=find(c>=s);
    x=x(1);
end
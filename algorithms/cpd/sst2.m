function [result]=sst2(x,w,n,l,g,m)
% Finds the points of change of x's dynamics. The output is returned as a
% row vector
% x = the input vector 
% w = the width of the dynamic window
% n = the number of windows to cosider in the past
% g = the width of the dynamic window in the future
% m = the number of windows to consider in the future

if(nargin<5)
    g=0;    
end
if nargin<6
    m=n;
end
   
echo off
H=zeros(w,n);
G=zeros(w,n);
c=numel(x);
result=zeros(1,c);
for t=w+n+1:1:c-w-m-1-g
    % construct the Herkel Matrix of the Past
    for nn=1:n
        H(:,nn)=x(t-w-nn+1:t-nn);
    end
    % calculate SVD of H
    [u,s,v]=svd(H,'econ');
    
    % find the matrix of the future G
    for nn=0:m-1
        G(:,nn+1)=x(t+g+nn:t+g+w-1+nn);
    end
    % find the eigen vectors of GG' in descending order
    G2=G*G';
    [v,d]=eig(G2);
    d2=diag(d);
    [dummy,di]=sort(d2,'descend');
    v=v(:,di);    
    u=u(:,1:l);
    beta=v(:,1);    
    result(t)=1-cos(subspace(u,beta)*2/pi);    
end

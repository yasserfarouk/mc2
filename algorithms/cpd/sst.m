% This file is a part of the MC2 toolbox developed by Y. Mohammand and T. Nishida.
%Please do not remove this comment
%
% Using this file is governed by the license of MC2 which you can find in LICENSE.md
% 
% You can find more information about this toolbox here:
% - Yasser Mohammad and Toyoaki Nishida, "MC2: An Integrated Toolbox for Change, Causality, 
%   and Motif Discovery", 29th International Conference on Industrial, Engineering & 
%   Other Applications of Applied Intelligent Systems (IEA/AIE) 2016, pp. 128 -- 141.
% - Yasser Mohammad and Toyoaki Nishida, "Data Mining for Social Robotics", Springer 2016.
%

function [result]=sst(x,w,n,l,g,m)
% Finds the points of change of x's dynamics. The output is returned as a
% row vector
%
% Inputs:
% =======
% x = the input vector 
% w = the width of the dynamic window
% n = the number of windows to cosider in the past
% g = the width of the dynamic window in the future
%
% Outputs
% =======
% m = the number of windows to consider in the future
%
% For more information please consult the following publications: 
% ===============================================================
% Tsuyoshi Ide and Keisuke Inoue Knowledge Discovery from Heterogeneous 
% Dynamic Systems using change point correlations, SIAM 2004, pp 571-575
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%

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
    alpha=(u'*beta);
    alpha=alpha./(norm(alpha,2));
    alpha2=zeros(size(beta));
    for i=1:l
        alpha2=alpha2+alpha(i).*u(:,i);
    end
    result(t)=1-alpha2'*beta;
%     if(result(t)>1)
%         result(t)=result(t)-floor(result(t));
%     end    
end

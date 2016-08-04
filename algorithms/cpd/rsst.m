function [result]=rsst(x,w,n,l,doPost,doThin)
% Finds the points of change of x's dynamics. The output is returned as a
% row vector
% x = the input vector 
% w = the width of the dynamic window
% n = the number of windows to cosider in the past
% l = the number of singular vectors to retain. if zero then RSST chooses
% doPost = if nonzero post processing to remove noise will be done
% doThin = if nonzero, a final thinning operation will be done
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, Robust Singular Spectrum Transform,
% The Twenty Second International Conference on Industrial, Engineering & 
% Other Applications of Applied Intelligent Systems (IEA/AIE 2009), 
% June 2009, Taiwan, pp 123-132
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%



percent=0.95;
if(nargin<3)
    n=w;
end
if(nargin<4)
    l=0;
end
if(nargin<5)
    doPost=0;
end
if(nargin<6)
    doThin=0;
end
m=n;
g=0;

noLFinding=l;
    
echo off
H=zeros(w,n);
G=zeros(w,n);
c=numel(x);
result=zeros(1,c);
noPast=0;
for t=w+n+1:1:c-w-m-1-g
    % construct the Hankel Matrix of the Past
    for nn=1:n
        H(:,nn)=x(t-w-nn+1:t-nn);
    end
    % calculate SVD of H
    [u,s,v]=svd(H,'econ');
    %find the number of eigenvectors to use
    if(noLFinding==0)
        if(s(1,1)<1e-9)
            u=u(:,1);
            l=1;
            noPast=1;
        else
            noPast=0;
            acc=cumsum(diag(s)); 
            acc=acc./max(acc);
            l=find(acc>percent,1,'first');
        end
    end
    % find the matrix of the future G
    for nn=0:m-1
        G(:,nn+1)=x(t+g+nn:t+g+w-1+nn);
    end
    % find the eigen vectors of GG' in descending order
    G2=G*G';
    [v,d]=eig(G2);
    d2=diag(d);
    [dd,di]=sort(d2,'descend');
    v=v(:,di);
    
    if(dd(1)<1e-9)
        result(t)=0;
    else
        if(noPast)
            result(t)=0;
        else
            if(noLFinding==0)
                acc=cumsum(dd); 
                acc=acc./max(acc);
                l2=find(acc>percent,1,'first');                
            else
                l2=l;
            end
            u=u(:,1:l);                        
            result(t)=0;
            for r=1:l2
                u2=v(:,r);
                result(t)=result(t)+1-cos(subspace(u,u2)*2/pi);
            end
            result(t)=result(t)/l2;
            
        end;
    end;
end
%do post processing using same w and n
if doPost~=0
    result=rsstpost(result,w,n);
end
%extra thining operation to keep only local minima
if doThin~=0
    result = thin(result);
end

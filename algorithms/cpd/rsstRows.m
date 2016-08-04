function [result]=rsst(x,w,n,l,doPost,doThin,justSST,g,m)
% Finds the points of change of x's dynamics. The output is returned as a
% row vector
% x = the input vector 
% w = the width of the dynamic window
% n = the number of windows to cosider in the past
% l = the number of singular vectors to retain. if zero then RSST chooses
% doPost = if nonzero post processing to remove noise will be done
% doThin = if nonzero, a final thinning operation will be done
% justSST = if nonzero then we will use RSST otherwise standard SST
% g = the width of the dynamic window in the future
% m = the number of windows to consider in the future

percent=0.95;
if(nargin<3)
    n=w;
end
if(nargin<4)
    l=0;
end
if(nargin<5)
    doPost=1;
end
if(nargin<6)
    doThin=1;
end
if(nargin<7)
    justSST=0;
end
if(justSST==0)
    m=n;
    g=0;
    l=1;
end
noLFinding=l;
    
echo off
H=zeros(n,w);
G=zeros(m,w);
c=numel(x);
result=zeros(1,c);
noPast=0;
for t=w+n+1:1:c-w-m-1-g
    % construct the Herkel Matrix of the Past
    for nn=1:n
        H(nn,:)=x(t-w-nn+1:t-nn);
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
        G(nn+1,:)=x(t+g+nn:t+g+w-1+nn);
    end
    % find the eigen vectors of GG' in descending order
    G2=G*G';
    [v,d]=eig(G2);
    d2=diag(d);
    [dd,di]=sort(d2,'descend');
    v=v(:,di);
    
    if(dd(1)<1e-9)
%         if(noPast)      %why am i doing that?
%             result(t)=1;
%         else
%             result(t)=0;
%         end;
        result(t)=0;
    else
        if(noPast)
            result(t)=0;
        else
            if(noLFinding==0)
                acc=cumsum(dd); 
                acc=acc./max(acc);
                l2=find(acc>percent,1,'first');
                %l=max(l,l2);
                %l2=l;
            else
                l2=l;
            end
            u=u(:,1:l);
            u2=v(:,1:l2);
            %u2=v2(1:l2,:)';
            result(t)=abs(subspace(u,u2)*2/pi);
            if(result(t)>1)
                result(t)=result(t)-floor(result(t));
            end
        end;
    end;
end
if doPost~=0
    result=rsstpost(result,w,n);
end
%extra thining operation to keep only local minima
if doThin~=0
%     r=result;
%     result2=diff(result);
%     result=zeros(numel(result2)+1,1);
%     result(2:numel(result))=result2;
%     result(result<0)=0;
%     result=result.*r;
    result = thin(result);
end

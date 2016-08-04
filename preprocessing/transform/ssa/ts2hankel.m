function H=ts2hankel(ts,L)
% converts a time series to a hankel matrix by embedding

if isvector(ts)    
    T=numel(ts);
    K=T-L+1;
    H=zeros(K,L);
    for k=1:K
        H(k,:)=ts(k:k+L-1);
    end
else
    T=size(ts,2);
    K=T-L+1;
    M=size(ts,1);
    H=zeros(K*M,L);
    for m=1:M
        for k=1:K
            H((m-1)*K+k,:)=ts(m,k:k+L-1);
        end
    end
end
end
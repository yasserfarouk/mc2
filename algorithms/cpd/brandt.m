function [c]=brandt(x,Nmax,L,p,q,step)
% finds change points in a 1D timeseries using brandt's method
    if ~exist('q','var')
        q=0;
    end
    if ~exist('step','var')
        step=1;
    end
    if isvector(x)
        x=x(:);
    end
    [T,~]=size(x);
    c=zeros(T,1);
    for t=min(Nmax,(L+p+q)):step:T
        N=min(Nmax,t);
        if q==0
            [~,~,~,~,~,sigma2All] = learnARLS(x(t-N+1:t,:),p);
            [~,~,~,~,~,sigma2Left] = learnARLS(x(t-N+1:t-L,:),p);
            [~,~,~,~,~,sigma2Right] = learnARLS(x(t-L+1:t,:),p);
        else
            [~,~,~,~,~,~,sigma2All] = learnARMALS(x(t-N+1:t,:),p,q);
            [~,~,~,~,~,~,sigma2Left] = learnARMALS(x(t-N+1:t-L,:),p,q);
            [~,~,~,~,~,~,sigma2Right] = learnARMALS(x(t-L+1:t,:),p,q);
        end
        sigma2All=sigma2All'*sigma2All;
        sigma2Left=sigma2Left'*sigma2Left;
        sigma2Right=sigma2Right'*sigma2Right;
        c(t)=-L*sqrt(sigma2Right)-(N-L)*sqrt(sigma2Left)+N*sqrt(sigma2All);
    end
end

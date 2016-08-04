function [g,Y,gDetails]=ccm(x,y,E, L,useHankelization)

% x,y if both are vectors input time-series to be tested for causality. 
%     if x is a matrix and y is [] then causality will be found between
%     columns of x
% E   embedding 
% Lmin,Lstep   [optional] the minimum value for time-series length used for
%           calculating ccm and its step up tp length(x)
if ~exist('useHankelization','var')
    useHankelization=false;
end
if ~isempty(y) && (~isvector(x) || ~isvector(y))
    error('either y should be empty or x and y should be vectors');
end
if ~isempty(y)
    x=[x(:),y(:)];
end

[T,N]=size(x);

if N<2
    error('at least 2 variables must be given');
end

if ~exist('L','var')
    L=T;
end
nTrials=floor(T/L);
xsave=x;
gDetails=zeros(N,N,nTrials);
if useHankelization 
    Y=zeros(N,N,L*nTrials);
else
    Y=zeros(N,N,(L-E+1)*nTrials);
end
epsilon=1e-6;
for t=1:nTrials
    x=xsave((t-1)*L+1:t*L,:);
    z=zeros(N,E,L-E+1);
    for i=E:L
        z(:,:,i-E+1)=x(i-E+1:i,:)';
    end
    
    for n=1:N
        xx=squeeze(z(n,:,:))';
        [idx,d]=knnsearch(xx,xx,'K',E+2);
        idx=idx(:,2:end); d=d(:,2:end);
        for i=1:L-E+1
            if d(i,end)<epsilon
                d(i,:)=ones(size(d(i,:)))./numel(d(i,:));
            elseif (d(i,1)/d(i,end))<epsilon
                d(i,:)=zeros(size(d(i,:)));
                d(i,1)=1;
            end
            d(i,:)=exp(-d(i,:)./(d(i,1))+epsilon);
            d(i,:)=d(i,:)/(sum(d(i,:))+epsilon);
        end
        for m=1:N
           if m==n
               continue;
           end
           if ~useHankelization
                y=zeros(1,L-E+1);
                for i=1:L-E+1
                    y(i)=d(i,:)*x(idx(i,:),m);
                end
                Y(m,n,(t-1)*L+E:t*L)=y;
           else
                yy=zeros(E,L-E+1);
                for i=1:L-E+1
                    yy(:,i)=squeeze(z(m,:,idx(i,:)))*d(i,:)';
                end
                y=hankel2ts(yy);
                Y(m,n,(t-1)*L+1:t*L)=y;
           end            
            if ~useHankelization
                gDetails(m,n,t)=corr(y(1,1:L-E+1)',x(E:L,m));
            else
                gDetails(m,n,t)=corr(y(1,1:L)',x(1:L,m));
            end
        end
    end
end
g=mean(gDetails,3);
end


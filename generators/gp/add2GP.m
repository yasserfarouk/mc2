function gp=add2GP(gp,X,y)
% Adds data to a GP
%
% % X is a D*n matrix where D is the number of dimensions of the input and n
%   is the number of example points (this is the form used by GP book by
%   Rasmussen and Wiliams)
n=size(X,2);
Kxoldxold=gp.Kxx;

Kxxold=gp.covFun(X,gp.X,gp.covFunParams);
Kxx=gp.covFun(X,X,gp.covFunParams);

Kxoldx=Kxxold';
Kxx=[Kxoldxold,Kxoldx;Kxxold,Kxx];
gp.X=[gp.X,X];
gp.y=[gp.y,y];
[u,s,v]=svd(Kxx);
s=diag(s);
logKxx=sum(log(s));
s=s.^-1;
Kxxinv=v*diag(s)*u';
n=size(X,2);
logpyX=-0.5*(gp.y*Kxxinv*gp.y'-logKxx-n*log(2*pi));
gp=struct('X',gp.X,'y',gp.y,'Kxx',Kxx,'Kxxinv',Kxxinv,'covFun',gp.covFun,...
    'covFunParams',gp.covFunParams,'sigmaN',gp.sigmaN,'logKxx',logKxx...
    ,'logpyX',logpyX,'fMean',gp.fMean,'fMeanX',[gp.fMeanX,gp.fMean(X)]);
end
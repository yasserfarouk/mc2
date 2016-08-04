function [fstar,m,K]=interpolateGP(Xstar,X,f,covFun,covFunParams,fMean)
% generates data at Xstar points given that we know the function values at X (i.e. f(X)=f)
%
% % X is a D*n matrix where D is the number of dimensions of the input and n
%   is the number of example points (this is the form used by GP book by
%   Rasmussen and Wiliams)
% covFun The covariance matrix to use. Default is covSE
% covFunParams The parameters of the covariance matrix to use in a cell
%               array containing the parameters of covFun in order. For
%               covSE (default covFun), this is {1} and sts l
% fMean         The mean of the GP functions (E(f(x))). Default is zero

if nargin<4
    covFun=[];
end
if nargin<5
    covFunParams=[];
end
if nargin<6
    fMean=0;
end

Kxsx=covFun(Xstar,X,covFunParams);
Kxx=covFun(X,X,covFunParams);
Kxxinv=inv(Kxx);
Kxsxs=covFun(Xstar,Xstar,covFunParams);
Kxxs=covFun(X,Xstar,covFunParams);
Ktmp=Kxsx*Kxxinv;
m=Ktmp*f;
K=Kxsxs-Ktmp*Kxxs;

fstar=mvnrnd(m,K)';

end
function [f,K]=generateGP(X,sigmaN,covFun,covFunParams,fMean)
% generates function values at columns of X given a covariance matrix from
% a GP
%
% X is a D*n matrix where D is the number of dimensions of the input and n
%   is the number of example points (this is the form used by GP book by
%   Rasmussen and Wiliams)
% sigmaN  The standard deviation of a Gaussian noise model to be added to the output
% covFun The covariance matrix to use. Default is covSE
% covFunParams The parameters of the covariance matrix to use in a cell
%               array containing the parameters of covFun in order. For
%               covSE (default covFun), this is {1} and sts l
% fMean         The mean of the GP functions (E(f(x))). Default is zero

if nargin<2 || isempty(sigmaN)
    sigmaN=0;
end

if nargin<3 || isempty(covFun)
    covFun=@covSE;
end

if nargin<4 || isempty(covFunParams)
    covFunParams={1};
end

D=size(X,1);
n=size(X,2);

if nargin<5 || isempty(fMean)
    fMeanX=zeros(1,n);
elseif isnumeric(fMean) && numel(fMean)==1
    fMeanX=fMean.*ones(1,n);
else
    fMeanX=fMean(X);
end
K=covFun(X,X,covFunParams);
K=K+sigmaN^2.*eye(size(K));
f=grand(fMeanX',K)';

end
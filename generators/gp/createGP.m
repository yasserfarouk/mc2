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

function gp=createGP(X,y,sigmaN,covFun,covFunParams,fMean)
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
if nargin<3 || isempty(sigmaN)
    sigmaN=0;
end
if nargin<4 || isempty(covFun)
    covFun=@covSE;
end
if nargin<5 || isempty(covFunParams)
    covFunParams=[];
end
n=size(X,2);
if nargin<6 || isempty(fMean)
    fMean=@meanLinear;
end

Kxx=covFun(X,X,covFunParams);
Kxx=Kxx+sigmaN^2.*eye(size(Kxx));
[u,s,v]=svd(Kxx);
s=diag(s);
logKxx=sum(log(s));
s=s.^-1;
Kxxinv=v*diag(s)*u';
n=size(X,2);
logpyX=-0.5*(y*Kxxinv*y'-logKxx-n*log(2*pi));
gp=struct('X',X,'y',y,'Kxx',Kxx,'Kxxinv',Kxxinv,'covFun',covFun,...
    'covFunParams',covFunParams,'sigmaN',sigmaN,'logKxx',logKxx...
    ,'logpyX',logpyX,'fMean',fMean,'fMeanX',fMean(X));

end
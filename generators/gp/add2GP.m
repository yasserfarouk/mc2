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
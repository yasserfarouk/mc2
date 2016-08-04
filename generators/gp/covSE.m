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

function k=covSE(Xp,Xq,hyp)
% generates a covariance matrix given a set of input vectors.
%
% X is a D*n matrix where D is the number of dimensions of the input and n
%   is the number of example points (this is the form used by GP book by
%   Rasmussen and Wiliams

% if nargin<2 || isempty(Xq)
%     Xq=Xp;
% end
if nargin<3 || isempty(hyp)
    hyp=[1,1];
end
if iscell(hyp)
    hyp=cell2mat(hyp);
end
l=hyp(1);
if numel(hyp)<2
    hyp(2)=1;
end
sigmaF=hyp(2);

%D=size(X,1);
n=size(Xp,2);
m=size(Xq,2);
k=eye(n,m);
l2=l^2;
for p=1:n;
    for q=1:m
        x=Xp(:,p)-Xq(:,q);
        k(p,q)=sigmaF.^2*exp(-0.5*sum(x.*x)/l2);
    end
end   
end
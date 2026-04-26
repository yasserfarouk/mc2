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

function show1DGaussianFun(m,Sigma,X,color,lightcolor)
% plots the mean with errorbars at the locations of X. Only works when X is a vector
    if ~isvector(X)
        error('X must be a vector');
    end
    n=numel(X);
    if n~=numel(m)
        error('m and X must be of the same length');
    end
    if size(Sigma,1)~=n || size(Sigma,2)~=n
        error('Sigma must be n*n if X is 1*n');
    end
    if nargin<4 || isempty(color)
        color=[0,0,0];
    end
    if nargin<5 || isempty(lightcolor)
        lightcolor = color + [0.6,0.6,0.6];
        lightcolor(find(lightcolor>1.0)) = 1.0;
    end
    
    hold on;
    
    X1=m(:)+2.*sqrt(diag(Sigma));
    X2=m(:)-2.*sqrt(diag(Sigma));
    X=X(:);
    C=[120,0,0];
     patch([X;X(end:-1:1)],[X1;X2(end:-1:1)], lightcolor, 'LineStyle', 'none');
     plot(X,m,'color',color,'LineWidth',3);
end
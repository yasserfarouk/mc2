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

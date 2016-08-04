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

function H=ts2hankel(ts,L)
% converts a time series to a hankel matrix by embedding

if isvector(ts)    
    T=numel(ts);
    K=T-L+1;
    H=zeros(K,L);
    for k=1:K
        H(k,:)=ts(k:k+L-1);
    end
else
    T=size(ts,2);
    K=T-L+1;
    M=size(ts,1);
    H=zeros(K*M,L);
    for m=1:M
        for k=1:K
            H((m-1)*K+k,:)=ts(m,k:k+L-1);
        end
    end
end
end
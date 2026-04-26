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

function [H,x]=hankelize(A,M)
% converts a mtrix to the corresponding hankel matrix by simple averaging
% M is the number of dimensions in the output time-series
if nargin<2
    M=1;
end
ITotal=size(A,1);
J=size(A,2);
I=ITotal/M;
T=I+J-1;
x=zeros(M,T);
for m=1:M
    H=A((m-1)*I+1:m*I,:);    
    x(m,1)=H(1,1);
    x(m,end)=H(end,end);
    for t=2:T-1
        if t<=J
            i=1; j=t;
        else
            i=t-J+1;j=J;
        end
        s=0;
        n=0;
        while(j>=1 && i<=I)
            s=s+H(i,j);
            n=n+1;
            j=j-1;
            i=i+1;
        end
        x(m,t)=s/n;
    end
end
H=ts2hankel(x,J);
end
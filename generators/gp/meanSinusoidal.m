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

function [ f ] = meanSinusoidal( X,hyp )
%MEANLINEAR Summary of this function goes here
%   Detailed explanation goes here
D=size(X,1);
n=size(X,2);
if nargin<1|| isempty(hyp)
    hyp=zeros(1+D*2,1);
end
% The hyberparameters are:
% hyp(1)=constant 
% hyp(2:D+1)=frequency
% hyp(D+2:end)=phase
f=hyp(2:D+1);
phase=hyp(D+2:end);
w0=hyp(1);
if numel(phase)<D
    phase=[phase,ones(1,D-numel(phase))];
end
    
    f=w0.*ones(1,n)+sin(2*pi.*f(:)'*X+phase(:)'*ones(D,n));

end


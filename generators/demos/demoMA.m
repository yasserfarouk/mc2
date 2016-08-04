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

figure; 
for i=1:8; 
    subplot(4,2,i); 
    a=[1,repmat(2*ones(1,3),1,i)]; 
    x=generateMA(a,1000,1); plot(x); 
    title(sprintf('$m=%d,\\left\\langle X\\right\\rangle=%5.3f,Var\\left(X\\right)=%5.3f$',numel(a),mean(x),var(x)),'interpreter','latex');  
end;
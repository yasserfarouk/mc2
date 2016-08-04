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

function [V,E,kind]=drawModel(M,n)
% draws a causality model using grPlot (M is the model and n (optional) the
% number of processes)
showDelay=0;
nSeries=size(M,1);
if(nargin<2); n=0; end;
n=max(nSeries,n);    
V=[cos(2*pi*[1:n]'/n) sin(2*pi*[1:n]'/n)]; 
kind='d';

E=[];
for i=1:nSeries    
    nI=size(M{i,1});
    for j=1:nI                
        if(showDelay)
            E=[E;M{i,1}(j,1),i,round(M{i,1}(j,2))];        
        else
            E=[E;M{i,1}(j,1),i];        
        end
    end
end
grPlot(V,E,kind,'%d','%d');
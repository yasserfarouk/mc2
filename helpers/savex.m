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

function savex(x,filename)
if nargin<2
    filename=sprintf('./Data/%s.dat',inputname(1));
end
f=fopen(filename,'w');
fprintf(f,'%s\n',inputname(1));
for i=1:numel(x)
    fprintf(f,'%8.5f\n',x(i));
end
fclose(f);
end
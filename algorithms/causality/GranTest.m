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

function [F1, F2]=grantest(u,y,m,n)
%GranTest Granger causality test

N = length(y);
if length(u)~=N, error('Improper input args.'), end
if nargin<4, n=m; end

thy_u = arx([y,u],[m, m, 0]);
thy = ar(y, m,'ls');
S1y=thy_u(1,1);
S2y=thy(1,1);

F1 = (N-2*m)*(S2y-S1y)/m/S1y;
Pf1=fpdf(F1,m,m);

thu_y = arx([u,y],[n, n, 0]);
thu = ar(u, n,'ls');
S1u=thu_y(1,1);
S2u=thu(1,1);

F2 = (N-2*n)*(S2u-S1u)/n/S1u;
Pf2=fpdf(F2,n,n);

if ~nargout
   if F2>F1
      disp('The second argument causes the first one.')
   else
      disp('The first argument causes the second one.')
   end
end

      
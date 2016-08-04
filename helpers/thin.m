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

function [result]=thin(data,th,varargin)
% applys a thinning operation on every row of its input
% data      the input row vector or array
% th        a threshold to stop operation should be a small value and must be
%           less than 1.0 (typical value is 1e-6)
if(nargin<2)
    th=1e-6;
end
n=size(data);
if(size(n)>2)
    error('The input must be either a single row vector or a two dimensional array');
end
% first find local maximums with some ridges
result=data;
for i=1:n(1);
    for j=2:n(2);
        if(result(i,j)>result(i,j-1));
            result(i,j-1)=0;        
        end;
    end
    for j=n(2):-1:2;
        if(result(i,j)<result(i,j-1));
            result(i,j)=0;        
        end;
    end
end
% changed=1;
% while(changed);
%     changed=0;
%     for i=1:n(1);
%         for j=1:n(2)-1;
%             if result(i,j)>th && (result(i,j)==result(i,j+1));
%                 result(i,j)=0;
%                 %result(i,j+1)=0;
%                 changed=1;
%             end;
%         end
%     end
% end;
% cmp=abs(result-data)>th;
% if max(cmp)>0.9;
%     result=thin(result,th);
% end
% second remove ridges to keep only local maximums

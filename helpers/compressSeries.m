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

function s=compressSeries(d,n,pointsAreRows)%,mx,mn)
% changes the series length to be equal to n. can be applied to vectors
% or 2d matrices
%
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%

if isvector(d)
    m=numel(d);
    if m ==n
        s=d;
    else
        if(size(d,1)==1)
            s=zeros(1,n);
        else
            s=zeros(n,1);
        end
        for i=1:n
            s(i)=d(max(1,floor(i*m/n)));
        end            
    end
else
    if nargin<3
        pointsAreRows=1;    
    end
    m=size(d,pointsAreRows+1);
    if m ==n
        s=d;
    else    
        if pointsAreRows
            s=zeros(size(d,1),n);
            for i=1:n
                s(:,i)=d(:,max(1,floor(i*m/n)));
            end    
        else
            s=zeros(n,size(d,2));
            for i=1:n
                s(i,:)=d(max(1,floor(i*m/n)),:);
            end    
        end
    end 
end
%s=((mx-mn)./(max(s)-min(s))).*s+mn;
end

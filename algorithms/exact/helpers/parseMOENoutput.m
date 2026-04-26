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

function [ locs,maxSmallDistance,means,stats ] = parseMOENoutput( str,x )
%Parses the output of any exact motif discovery algorithm in the toolbox
%   
    A=str2num(str);    
    n=size(A,1);
    locs=cell(n,1);
    if(nargout>2)
        means=cell(n,1);
        stats=cell(n,1);
    end
    maxSmallDistance=max(A(:,4));
    A=A(:,1:3);
    A=int32(A);
    for i=1:n
        locs{i}=[A(i,2),A(i,1)+A(i,2)-1;...
                 A(i,3),A(i,1)+A(i,3)-1];
       if(nargout>1)
            [stats{i},means{i}]=calcStatsAndMean(locs{i},x,@abs);
       end
    end   
    
end


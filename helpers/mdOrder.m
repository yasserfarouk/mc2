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

function  [locs,means,stats]= mdOrder(locs,means,stats,order)
% orders motifs by the length*nOccurrences then first occurrence of each of them and orders
% occurrences witin each motif by their occurrence 
%
% if order is negative it is sorting in descending order otherwise
% assending order. Default is ascending

    if nargin<nargout
        error('number of inputs must be larger than or equal the number of outputs');
    end
    if nargin<4
        order=1;
    else
        if order==0
            order=1;
        end
    end
    order=sign(order);
    nLocs=length(locs);    
    strt=zeros(nLocs,2);
    for i=nLocs:-1:1
        locs{i}=sortrows(locs{i});       
        strt(i,1)=(locs{i}(1,2)-locs{i}(1,1)).*(size(locs{i},1)-1);
        strt(i,2)=locs{i}(1,1);
    end
    [dummy,I]=sortrows(strt,order);
    locs=locs(I);
    if nargout>1
        means=means(I);
    end
    if nargout>2
        stats=stats(I);
    end
end

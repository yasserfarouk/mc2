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

function stop=isMeanDistBetweenOccurencesSmall(x,newbegends,maxSmallDist...
    ,minLargeDist,dFun,dFunParams)
% Tests the smallness of a distance. (internal function)
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%

    n=size(newbegends,1);
    meanDist=[]; l=0; 
    for i=1:n
        for j=i+1:n
            dist=dFun(x(newbegends(i,1):newbegends(i,2)),...
                      x(newbegends(j,1):newbegends(j,2)),dFunParams);
            if isempty(meanDist)
                meanDist=dist; maxDist=dist;
                l=1;
            else                    
                maxDist=max([maxDist,dist]);
                if maxDist>minLargeDist
                    stop=1;
                    return;
                end
                meanDist=(l*meanDist+dist)/(l+1);                    
                l=l+1;                    
            end
        end
    end
    stop=(meanDist>maxSmallDist);        
end
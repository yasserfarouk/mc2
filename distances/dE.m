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

function d=dE(s1,s2,optional)
% Distance with possible normlaization. 
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
    normalize=0;
    power=2;
    ignoreDC=1;
    minRatio=0.001;

    if nargin>2 && ~isempty(optional)
        normalize=optional{1};
        if length(optional)>1
            power=optional{2};
        end
        if length(optional)>2
            ignoreDC=optional{3};
        end
        if length(optional)>3
            minRatio=optional{4};
        end
    end    
    switch (normalize) 
        case 1 % normalize by subtracting the mean
            s1=s1-mean(s1).*ones(size(s1));
            s2=s2-mean(s2).*ones(size(s2));
        case 2  % normlaize by range
            delta=max(s1)-min(s1);
            if delta>minRatio
                s1=(s1-mean(s1).*ones(size(s1)))./(delta);
            else
                if ignoreDC
                    d=inf;
                    return;
                end
            end        
            delta=max(s2)-min(s2);
            if delta>minRatio
                s2=(s2-mean(s2).*ones(size(s2)))./(delta);
            else
                if ignoreDC                
                    d=inf;
                    return;
                end
            end               
        case 3 % normalize by mean and stadard deviation
            s1=zscore(s1);
            s2=zscore(s2);            
        otherwise     
    end
    d=mean(abs((s2(:)-s1(:)).^power));
end

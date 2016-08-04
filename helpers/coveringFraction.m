function [fraction,reverse,coveredRange]=coveringFraction(covered,covering)
%finds how much of covered is covered by covering (internal function)
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
    if covered(1)>covering(2) || covering(1)>covered(2)
        fraction =0;
        reverse=0;
        coveredRange=[];
        return;
    end
    if (covering(1)<=covered(1) && covering(2)>= covered(2))
        fraction =1;        
        reverse=(covered(2)-covered(1))./(covering(2)-covering(1));        
        coveredRange=covered;
        return;
    end
    coveredRange=[max([covered(1),covering(1)]),min([covered(2),covering(2)])];
    fraction=(coveredRange(2)-coveredRange(1)+1)/(covered(2)-covered(1)+1);
    reverse=(coveredRange(2)-coveredRange(1)+1)/(covering(2)-covering(1)+1);    
end
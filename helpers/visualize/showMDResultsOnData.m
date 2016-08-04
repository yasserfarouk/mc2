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

function showMDResultsOnData(x,locs,useMotifRaising)
if ~exist('useMotifRaising','var') || isempty(useMotifRaising)
    useMotifRaising=true;
end
[T,m]=size(x);
plot(x,'g','LineWidth',0.5);
hold on;
n=numel(locs);
names=cell(n,1);
if useMotifRaising
    smallValue=(max(x)-min(x))*0.1;
    raisedValues1=smallValue.*(0:(max([2,ceil(n/2)])));
    raisedValues2=-raisedValues1(2:end);
    raisedValues=zeros(numel(raisedValues1)+numel(raisedValues2),1);
    raisedValues(2:2:end)=raisedValues1(2:end);
    raisedValues(3:2:end)=raisedValues2;
else
    raisedValues=zeros(n,1);
end
for i=1:n
    begend=locs{i};
    names{i}=sprintf('Motif %d',i);
    y=nan(T,m);
    K=size(begend,1);
    for k=1:K
        y(begend(k,1):begend(k,2))=x(begend(k,1):begend(k,2))+raisedValues(i);
    end
    plot(y,'LineWidth',2);
end


end
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

function th=selectThreshold(x,L,distFun)
if nargin<3 || isempty(distFun)
    distFun=@dE;
end
T=numel(x);

nSamples=T;
fraction=1/30.0;
minSeparation=L;
nTrials=3;
tmpLoc=randi(T-L+1,nSamples,1);
tmpLoc2=randi(T-L+1,nSamples,1);
indx=find(abs(tmpLoc-tmpLoc2)<minSeparation);
for i=1:nTrials
    
    if isempty(indx)
        break;
    end
    nBadSamples=numel(indx);
    tmpLoc(indx)=randi(T-L+1,nBadSamples,1);
    tmpLoc2(indx)=randi(T-L+1,nBadSamples,1);
    indx=find(abs(tmpLoc-tmpLoc2)<minSeparation);
end

tmpLoc(indx)=[];
tmpLoc2(indx)=[];
nSamples=numel(tmpLoc);
d=zeros(nSamples);
for i=1:nSamples
    d(i)=dE(x(tmpLoc(i):tmpLoc(i)+L-1),x(tmpLoc2(i):tmpLoc2(i)+L-1));
end
d=sort(d);
th=d(ceil(fraction*T));
end
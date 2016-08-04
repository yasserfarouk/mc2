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

function [ x ] = generateRandomWalk(T,n,delta,p,x0,useGaussian )
%Generate data from a random walk. Timeseries are returned in column vectors
%   T       Time series length
%   n       number of time-series
%   delta   if >0 amount of increase/decrease at evey step
%           if <0 the maximum amount of increase/decrease at every step
%           multiplied by -1. The actual incrase/decrease at every step is
%           choosen randomly from a uniform distribution with a minimum of
%           zero and this maximum
%   p       probability of increasing
%   x0      initial value
%   useGaussian   if true a standard Normal distribution will be used for
%           deciding whether to increase or decrease the time-series value
% output
% ======
%   x       T*n matrix containing n time-series of length T
if ~exist('n','var') || isempty(n)
    n=1;
end
if ~exist('delta','var') || isempty(delta)
    delta=1;
end
if ~exist('p','var') || isempty(p)
    p=.5;
end
if ~exist('x0','var') || isempty(x0)
    x0=0;
end
if ~exist('useGaussian','var') || isempty(useGaussian)
    useGaussian=false;
end
useRandStep=false;
if delta<0
    delta=-delta;
    useRandStep=true;
end
x=zeros(T,n);
x(1,:)=x0;
if useGaussian
    pt=randn(T,n);
else
    pt=rand(T,n);
end
if useRandStep
    d=rand(T,n).*delta;
else
    d=delta.*ones(T,n);
end
for t=2:T    
    x(t,:)=x(t-1,:)+((-1).^(p-pt(t,:)<0)).*d(t,:);    
end

end


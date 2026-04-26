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

function x=generateMarkovChain(mu0,sigma0,sigma,T,a)
% generates data from a markov chain using p0 as initial probability and A
% as the transition matrix
%
% mu0       An N*1 vector giving the mean of the initial distribution used
%           for the first value of x
% sigma0       An N*N matrix giving the covariance of the initial distribution used
%           for the first value of x
% sigma     An N*N matrix giving the covariance of the initial distribution used
%           for x_t
% T         The lenght of the output time-series
% a     [optional] The parameters of the MA(m) process and m=numel(a) in the order
%       [a_0, a_1, ... a_m]
N=numel(mu0);
if ~exist('a','var') || isempty(a)
    a=[];
end
addNoise=~isempty(a);
assert(all(size(sigma0)==[N,N]));
assert(all(size(sigma)==[N,N]));
if addNoise
    ma=generateMA(a,T,N);
else
    ma=zeros(T,N);
end
if N==1
    x=zeros(T,1);
    r=randn(T,1);
    x(1)=mu0+sigma0*r(1)+ma(1);
    
    for t=2:T        
        x(t)=(x(t-1))+sigma*r(t)+ma(t);
    end
else
    x=zeros(T,N);    
    x(1,:)=mvnrnd(mu0',sigma0)+ma(1,:);
    for t=2:T        
        x(t,:)=mvnrnd(x(t-1,:)',sigma)+ma(t,:);
    end
end

end

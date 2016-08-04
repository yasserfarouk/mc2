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

function x=generateAffineStateSpace(A,B,C,sigma_epsilon, sigma_lambda,T,s0)
% generates data from an affine state space model
%
% The model used for generation has the form:
%
% s_t=As_{t-1}+B*epsilon
% x_t=Cs_t+lambda
%
% where epsilon and lambda are generated from zero-mean Gaussians with the
% given covariance matrices. If the covariance matrices are not given they
% are assumed to equal I.
%
%
Ns=size(A,1);
Nx=size(C,1);
if ~exist(sigma_epsilon,'var') || isempty(sigma_epsilon)
    sigma_epsilon=eye(Ns,Ns);
end
if ~exist(sigma_lambda,'var') || isempty(sigma_lambda)
    sigma_lambda=eye(Nx,Nx);
end
if ~exist(s0,'var') || isempty(s0)
    s0=zeros(Ns,1);
end
Neps=size(sigma_epsilon,1);
Nlam=size(sigma_lambda,1);
assert(size(B,1)==Ns & size(B,2)== Neps,'dimensions of B are not correct');
assert(size(C,2)== Nlam,'dimensions of B are not correct');
x=zreos(T,Nx);
for t=1:T
    s0=A*s0+B*mvnrnd(zeros(Neps,1),simga_epsilon);
    x(t,:)=C*s0+mvnrnd(zeros(Neps,1),simga_lambda);
end

end
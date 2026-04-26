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

function [x,s,y,Sigma_y]=generateGMR(prior,mu,sigma,T,a)
% generates data from a GMR process
%    
% p0        an K elements vector representing the priors of K gaussians. It must sum to 1 
% mu        An N+1*Ks matrix where mu(:,k) gives the mean of the kth Gaussian
%           (N data dimensions plus time)
% sigma     An N+1*N+1*Ks matrix giving the covariances of the K Gaussians
% T         The lenght of the output time-series
% a     [optional] The parameters of the MA(m) process and m=numel(a) in the order
%       [a_0, a_1, ... a_m]
if ~exist('a','var') || isempty(a)
    a=[];
end
n=size(mu,1)-1;
K=size(mu,2);
if ~exist('prior','var') || isempty(prior)
    prior=ones(K,1)/K;
end
if ~exist('a','var')
    a=[];
end
s=zeros(T,n);
if isempty(a)
    ma=zeros(T,n);
else
    ma=generateMA(a,T,n);
end

tt=1:T;
Pxi=zeros(T,K);
for i=1:K
  Pxi(:,i) = prior(i).*gaussPDF(tt, mu(end,i), sigma(end,end,i));
end
beta = (Pxi./repmat(sum(Pxi,2)+realmin,1,K))';
%% Compute expected output distribution, given input tt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = zeros(n, T);
Sigma_y = zeros(n, n, T);
for i=1:T
  % Compute expected means y, given input tt
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for j=1:K
    yj_tmp = mu(1:n,j) + (sigma(1:n,end,j)/(sigma(end,end,j))) * (tt(i)-mu(end,j));
    y(:,i) = y(:,i) + beta(j,i).*yj_tmp;
  end
  % Compute expected covariance matrices Sigma_y, given input x
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for j=1:K
    Sigmaj_y_tmp = sigma(1:n,1:n,j) - ((sigma(1:n,end,j)/(sigma(end,end,j))))*sigma(end,1:n,j);
    Sigma_y(:,:,i) = Sigma_y(:,:,i) + beta(j,i)^2.* Sigmaj_y_tmp;
  end
  s(i,:)=mvnrnd(y(:,i),Sigma_y(:,:,i))';
end
x=s+ma;
y=y';
end


function prob = gaussPDF(Data, Mu, Sigma)
%
% This function computes the Probability Density Function (PDF) of a
% multivariate Gaussian represented by means and covariance matrix.
%
% Author:	Sylvain Calinon, 2009
%			http://programming-by-demonstration.org
%
% Inputs -----------------------------------------------------------------
%   o Data:  D x N array representing N datapoints of D dimensions.
%   o Mu:    D x K array representing the centers of the K GMM components.
%   o Sigma: D x D x K array representing the covariance matrices of the 
%            K GMM components.
% Outputs ----------------------------------------------------------------
%   o prob:  1 x N array representing the probabilities for the 
%            N datapoints.     

[nbVar,nbData] = size(Data);

Data = Data' - repmat(Mu',nbData,1);
prob = sum((Data*inv(Sigma)).*Data, 2);
prob = exp(-0.5*prob) / sqrt((2*pi)^nbVar * (abs(det(Sigma))+realmin));

end
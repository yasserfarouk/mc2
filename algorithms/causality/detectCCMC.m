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

function [M,G]=detectCCMC(x,alpha,E,maxSteps,maxT,varargin)
% applies CCM causality test to all the dimensions of X and generates the causality model
% from which this data was generated
%
%
% Inputs:
% ======
% x             n*T array of time series to find the relations between its
%               n timeseries
% alpha         The significance level to be used to decide causality
%               (number between zero and one but should be near one)
% E             Embedding dimensions
% maxSteps      maximum number of lengths to use
% maxT          maximum length to use 
% varargin      see the switch statement
% 
% Outputs:
% =======
% M             The model (look at produceMulti to understand its
%               structure)

if nargin<2
    alpha=0.8;
end
if nargin<3
    E=2;
end
if nargin<4
    maxSteps=100;
end
[n,T]=size(x);
if nargin<5
    maxT=T;
end
step=max(E,floor((maxT-E+1)/maxSteps));

L=[E*2:step:maxT];
if L(end)~=maxT
    L(end)=maxT;
end
nTrials=numel(L);
useHankelization=true;

G=zeros(n,n);

rho=zeros(n,n,nTrials);

for t=1:nTrials
    g=ccm(x',[],E,L(t),useHankelization);
    rho(:,:,t)=g;
end
for i=1:n
    for j=1:n
        G(i,j)=rho(i,j,end);
        if G(i,j)<alpha
            G(i,j)=0;
        end
    end
end

M=cell(n,5);
for i=1:n
    parents=find(G(n,:)>0.0001);
    np=numel(parents);
    M{i,1}=cell(np,5);
    for k=1:np
        M{i,1}{k,1}=parents(k);
        M{i,1}{k,2}=nan;
        M{i,1}{k,3}=nan;
        M{i,1}{k,4}=G(i,parents(k));
        M{i,1}{k,5}=nan;
    end
    M{i,2}=0;
    M{i,3}=[];
    M{i,4}=1;
end
end
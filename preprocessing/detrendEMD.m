function [y,Tr]=detrendEMD(x,t)
% detrends the time-series x using EMD
%
% x     A vector representing a time-series
% t     The independent variable. Must be a vector of the same length as x

% y     The time-sereis after detrending
% Tr    The estimated Trend

assert(isvector(x),'x must be a vector');

T=numel(x);
if ~exist('t','var') ||  isempty(t)
    t=1:T;
end

IMF=emd(x);
Tr=IMF(end,:);
if size(x,1)~=1
    Tr=Tr';
end
y=x-Tr;
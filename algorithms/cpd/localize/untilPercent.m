function i=untilPercent(x,varargin)
% returns the index at which the percent of power specified as the first
% optional parameter is achieved.  default is 95
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
nArgs=size(varargin,2);
if(nArgs<1)
    p=0.95;
else
    if(varargin{1}>1)
        p=varargin{1}/100;
    else
        p=varargin{1};
    end
end
n=numel(x);
m=sum(x)*p;
s=0;
for i=1:n
    s=s+x(i);
    if abs(s-m)<eps
        return;
    elseif s>m
        return;
    end
end
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

function y=rsstpost(x,w,n,C,em,es)
% applies rsst post processing
%
% x         output of RSST
% w         the window width
% n         the number of windows used in RSST
% C         constant multiplied by the difference in trend (default 1)
% em        exponent of mean effect (default 1.5)
% es        exponenet of standard deviation difference (default 0)
% y         output of this function after post processing
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, Robust Singular Spectrum Transform,
% The Twenty Second International Conference on Industrial, Engineering & 
% Other Applications of Applied Intelligent Systems (IEA/AIE 2009), 
% June 2009, Taiwan, pp 123-132
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%
if(nargin<4)
    C=1;
end
if(nargin<5)
    em=1.5;
end
if(nargin<6)
    es=0;
end
l=numel(x);
y=zeros(size(x));
width=round(w/2);
for i=max(w+n+1,width):l-width-1
    s=x(i+1:i+width);
    mf=mean(s);
    segf=std(s);
    p=x(i-width:i-1);
    mp=mean(p);
    segp=std(p);
    y(i)=x(i)* C * ((abs(mf-mp))^em) *((abs(segf-segp))^es);
end
y=y./max(y);


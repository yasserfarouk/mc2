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

function dist=KLDiv(P,Q)
%  dist = KLDiv(P,Q) Kullback-Leibler divergence of two discrete probability
%  distributions
%  P and Q  must be vectors (either column or row but same)

if size(P,2)~=size(Q,2) || (size(P,1)~=size(Q,1))
    error('the number of columns in P and Q should be the same');
end
if sum(~isfinite(P(:))) + sum(~isfinite(Q(:)))
   error('the inputs contain non-finite values!') 
end
epsilon=1e-9;
% add small epsilon to get rid of zeros
if(numel(find(Q==0))~=0)
    Q=Q./max(Q);
    Q=Q+epsilon;
end

% normalizing the P and Q
Q = Q ./sum(Q);
P = P ./sum(P);
temp =  P.*log(P./Q);
temp(P==0)=0;
%temp(isnan(temp))=0;% resolving the case when P(i)==0
dist = sum(temp);   
end



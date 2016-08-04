function q=aesr(x,t,nPoints)
% finds the Accomulative Equivalence Sampling Rate between x and t pdfs given that t
% represents true probability of occurrence of some event.
%
% ESR is the probability that if I sample from x the sample will be within
% nPoints from some nonzero sample in t. The accomuative version
% accomulates these results over all nPoints<=the given
%
% x         output of change point algorithm. If x is n*m matrix then it is
%           treated as m outputs of m algorithms with a time series of
%           length n and nPoints must be n points. This means the each
%           algorithm's output must be in a column not a row
% t        true change points (e.g. 1 at change points and 0 otherwise)
%           the value at each change point may be treated as weight and the
%           more this value is the more important this change point will be
%           in the result. This must be a vector
% nPoints   width around change points that is assumed OK. 
%
%
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, On Comparing SSA-based Change Point
% Discovery Algorithms, IEEE SII 2011, 938-945 
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%


q=0;
for i=1:nPoints+1
    q=q+esr(x,t,i-1);
end
q=q./(nPoints+1);
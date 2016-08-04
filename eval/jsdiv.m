function r=jsdiv(p,q)
% computes Jennsen-Shannon Divergence but with correction in case both pi
% and qi are zeros
%


t=2./(p+q);
r=p*log2(p.*t)+q*log2(q.*t);
end
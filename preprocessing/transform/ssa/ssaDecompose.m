function [h,zeta,U,S,V]=ssaDecompose(ts,L,minZeta)
% Applies SSA decomposition returning a cell array of decomposition
% matrices and corresponding weights
if nargin<3
    minZeta=1e-6;
end
% if isvector(ts)
%     H=ts2hankel(ts,L);
% else    
%     M=size(ts,1);
%     H=[];
%     for m=1:M
%         Hm=ts2hankel(ts,L);
%     end
%     H=[H,Hm];
%     
% end
H=ts2hankel(ts,L);
[U,s,V]=svd(H,'econ');
S=diag(s);
S(S<minZeta)=[];
d=numel(S);
h=cell(d,1);
zeta=S./sum(S);
for i=1:d
    h{i}=S(i).*(U(:,i)*V(:,i)');
end
end  
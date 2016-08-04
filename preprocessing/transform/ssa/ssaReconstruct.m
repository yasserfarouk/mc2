function [x,zetag,Hgs]=ssaReconstruct(h,zeta,groups,nDimensions)
% reconstructs the expansion time-series from a cell array of expansion
% matrices and their corresponding zeta and user defined grouping
%
%
% h         cell array of matrices found by ssaDecompose
% zeta      expansion matrix weights found by ssaDecompose
% groups    One of the following:
%            Option 1: nG*1 cell array each contains a vector of integers representing
%                       the members of a group. The user MUST be sure that the groups are mutually exclusive
%            Option 2: A vector of integers containing the number of
%                       expansion matrices to be combined in each group in
%                       order
%            Option 3: A float giving the total weight to put in the first
%            group. The rest is put in a second group. This leads to two
%            groups always
%            Option 4: nothing or null is equivalent to option 3 with 0.99
% nDimensions The number of dimensions in the output time-series. 1 if
%           omitted
% Output
% ======
% x         A cell array containing nG time-series comprising the expansion
% zeta      The weight of each expansion time-series

if nargin<3
groups=0.99;
end
if nargin<4
    nDimensions=1;
end
nH=numel(h);

if ~iscell(groups)
    
    if not (isvector(groups) && numel(groups)>1)
        sz=cumsum(zeta);
        g=[nH,0];
        for i=1:numel(sz)
            if sz(i)>=groups
                g=[i,nH-i];
                break;
            end;
        end
        groups=g;
    end
    if ~isvector(groups)
        error('invalid input for groups');
    end
    nG=numel(groups);
    g=cell(nG,1);
    start=1;
    for i=1:nG
        g{i}=start:start+groups(i)-1;
        start=g{i}(end)+1;
    end
    groups=g;
end

nG=numel(groups);



zetag=zeros(nG,1);
Hgs=cell(nG,1);
x=cell(nG,1);
for g=1:nG
    Hgs{g}=zeros(size(h{1}));
    for i=groups{g}
        Hgs{g}=Hgs{g}+h{i};
        zetag=zetag+zeta(i);
    end
end

for g=1:nG
    x{g}=hankel2ts(Hgs{g},nDimensions);
end

end
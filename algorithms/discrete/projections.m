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

function [locs,means,stats] = ...
    projections(xD,nIterations,fractionKept,x,nPointsPerSymbol,R,varargin)
% discovers motifs in discretized time series using the projections
% algorithm
%
% INPUTS:
% =======
% xD        The discretized time series as M*L matrix where M is the number
%           of strings in the discretized version and L is the length of
%           each string. It assumes that the input string is obtained using
%           SAX or a similar algorithm from a continuous signal (x) and
%           that xD represents all substrings of length L of the
%           discretized version. in SAX this can be obtained by calling:
%           xD=timeseries2symbol(x,nPointsPerSymbol*L,L,alphabit_size).
%           Notice that the length of the motif to be found is L in
%           discrete space or L*nPointsPerSymbol in real space
% nIterations   The number of iterations of projections to be carried out
% fractionKept  The fraction of the maximum number of collisions to be
%               considered important and kept in the collision matrix
% x         the original time series before disretization. Will
%           be used only if outputs from locs forward are needed
% nPointsPerSymbol  The number of points represented by each symbol
% R        The maximum allowable distance between any two occurrences of a
%           motif per point.
% varargin  see the switch statement
%
% OUTPUTS:
% ========
% locsD     a m*1 cell array representing the row-indices of m motifs in the
%           xD array. each element of the cell-array contains a n_i*1
%           array representing the row representing the occurrence
% motifsD   The means of the motifs found as a m*1 cell array each with a
%           mean of one motif in the discretized space
%
% locs      a m*1 cell array representing the locations of m motifs in the
%           time series. each element of the cell-array contains a n_i*2
%           array representing the beginning and end of one occurrence of
%           the motif
% means     The means of the motifs found as a m*1 cell array each with a
%           mean of one motif
% stats     the statistics of each motif in the form of a structure array.
%           members:
%               mean    The mean difference between corresponding
%                           points in the occurrences
%               var     The standard deviation difference between corresponding
%                           points in the occurrences
%               max     The maximum difference between corresponding
%                           points in the occurrences
%               min     The minimum difference between corresponding
%                           points in the occurrences
%               median  The mmedian difference between corresponding
%                           points in the occurrences
%               mode    The mode of the difference between corresponding
%                           points in the occurrences
%
% Notes:
% ======
% - Inputs x and PointsPerSymbol must be supplied if outputs locs, means,
%   stats (or any of them) are needed
%

[n,strLen]=size(xD);

visualize=0;
localizationMethod='chiu';
nMinHits=1;             % used only with localizationMethod chiu.
                        % the minimum number of hits within the range
                        % for a motif
nSymbolsInHash=max(1,round(strLen*.7));
dFun=@dE;
dPointDiffFun=@abs;
removeOverlappingOccurrences=0;
% initialize the sparse matrix M to have no values at all (all zeros)
%M=sparse(n,n);
Mr=[]; Mc=[]; Mv=[];

% we fill only the lower triangle of the collision martix because we know
% that it is symmetric. We also do not fill the diagonal because it
% represents the trivial match between a pattern and itself
for iteration=1:nIterations
    % for each iteration we select a random set of points in the string
    indx=randi(strLen,nSymbolsInHash);
    % we find the hash of the string as its sum
    hashs=sum(xD(:,indx),2);
    % makes sure that minimum hash is 1 to reduce space
    hashs=hashs-min(hashs)+1;   
    % find for each pssible hash value all the rows that hashed to it
    nBins=max(hashs);
    bins=cell(nBins,1);    
    for i=1:n
        bins{hashs(i),1}=[bins{hashs(i),1},i];
    end
    % we now find collisions defined as rows that hashed to the same
    % location. Notice that c will always be less than r and this is how we
    % are filling only the lower triangle under the diagonal
    for i=1:nBins
        nColliding=numel(bins{i,1});
        for c=1:nColliding-1
            cc=bins{i,1}(c);          % cc, bb are the rows that hashed to i
            for b=c+1:nColliding
                bb=bins{i,1}(b);
                if abs(bb-cc)>=strLen % ignore overlapping windows
                    %M(bb,cc)=M(bb,cc)+1;                    
                    % add one collision to the sparse M matrix
                    Mr=[Mr,bb];
                    Mc=[Mc,cc];
                    Mv=[Mv,1];
                end
            end
        end
    end
end

M=sparse(Mr,Mc,Mv,n,n);
% find the maximum collisions entry with its location (r,c)
[mxC,cs]=max(M,[],2);
[mx,r]=max(mxC);
c=cs(r);
% remove elements of the matrix that are less than the specified threshold
if visualize
    figure; 
    subplot(1,2,1);
    colormap(gray); 
    image( -M,'CDataMapping','scaled');
end
M(M<mx*fractionKept)=0;
if visualize
    subplot(1,2,2); 
    colormap(gray); 
    image( -M,'CDataMapping','scaled');
end
[Mr,Mc,Mv]=find(M);
nM=numel(Mr);

if exist('R','var')
    % do check with original time series data
    if ~isempty(R)
        for i=1:nM
            first=x(Mr(i):Mr(i)+nPointsPerSymbol*strLen-1);
            second=x(Mc(i):Mc(i)+nPointsPerSymbol*strLen-1);
            d=dFun(first,second);
            if d>R
                Mr(i)=0;
                Mc(i)=0;
                %Mv(i)=0;            
            end
        end    
        Mr(Mr==0)=[];
        Mc(Mc==0)=[];
        nM=numel(Mr);
    else
        % TODO implement minnen's method here to find R then use it as
        % above
    end
    %Mv(Mv==0)=[];
end
% find the subsequences with maximum number of near subsequences
switch( localizationMethod)
case 'chiu'
    % find the top subsequences in terms of the number of similiar
    % subsequences. Remove anything that has less than nMinHits hits.
    % Finally generate the output locations of motifs.
    hits=zeros(n,1);
    hitLocs=cell(n,1);
    for i=1:nM
        hits(Mr(i))=hits(Mr(i))+1;
        hitLocs{Mr(i)}=[hitLocs{Mr(i)},Mc(i)];
        hits(Mc(i))=hits(Mc(i))+1;
        hitLocs{Mc(i)}=[hitLocs{Mc(i)},Mr(i)];
    end
    [hits,motifLocs]=sort(hits);
    motifLocs(hits<nMinHits)=[];
    nM=numel(motifLocs);
    locs=cell(nM,1);    
    for i=1:nM        
        myhitLocs=hitLocs{motifLocs(i)};
        nHits=numel(myhitLocs);
        locs{i}=zeros(nHits+1,2);
        locs{i}(1,:)=[motifLocs(i),motifLocs(i)+nPointsPerSymbol*strLen-1];   
        for j=1:nHits
            if ~removeOverlappingOccurrences || abs(myhitLocs(j)-locs{i}(j,1))>=strLen
                locs{i}(j+1,:)=[myhitLocs(j),myhitLocs(j)+nPointsPerSymbol*strLen-1];            
            end
        end
        locs{i}(locs{i}(:,1)==0,:)=[];
    end
case 'all'
    locs=cell(nM,1);    
    for i=1:nM
        locs{i}=[Mc(i),Mc(i)+nPointsPerSymbol*strLen-1;Mr(i), Mr(i)+nPointsPerSymbol*strLen-1];   
    end
end

if nargout==2
    means=cell(nM,1);
end
if nargout>2
    means=[];
    stats=[];
end
for i=1:nM
    if nargout>2
        [stats,means]=updateStatsAndMean(locs{i},x,dPointDiffFun,stats,means,i);
    else
        if nargout>1
            first=x(locs{i}(1,1):locs{i}(1,2));
            second=x(locs{i}(2,1):locs{i}(2,2));
            means{i}=mean([first(:),second(:)],2)';
        end
    end
end

end

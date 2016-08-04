function [coveringNone,coveringPartial,coveringMulti,covered,extras,badMotif,covering,coveredby]...
    =mdquality(locs,locsT,minCover2Consider,maxOverlap2Neglect)
% calculates the motif discovery quality statistics
%
%
% INPUTS:
% =======
% locs              Discovered motif locations in a cell array as outputed 
%                   from dgr()
% locsT             ground truth location of motifs in the same format as 
%                   locs
% minCover2Consider [optiona] the minimum percent of coverage to be 
%                    considered  (default 0)
% 
% maxOverlap2Neglect [Optional] maximum coverage to be considered when
%                    calculating the found motifs that cover multiple ground
%                    truth motifs (default 0)
%
% OUTPUTS:
% =======
% coveringNone      The fraction of the found motif that cover nothing in
%                   the ground truth
% coveringPartial   The fraction of the found motif where only some of the
%                   occurrences cover ground-truth occurrences but the rest
%                   is covering nothing
% coveringMulti     The fraction of the found motif that cover multiple
%                   motifs in the groung truth
% covered           The fraction of each ground truth motif that is covered
%                   by found motifs corresponding to this motif
% extras            The length of all extra parts in found covering
%                   occurrences for a motif as a fraction of this motif's
%                   total occurrence length
% badMotif          a vector with nonzero for each found motif that is
%                   either not covering anything (1), covering partially (2)
%                   or covering multiple ground truth motifs (3)
% covering          a cell array with an element for each motif in the
%                   found motifs. this element is an array of 7 columns
%                   giving the following information in order:
%                   1: the occurrence of the found motif
%                   2: the ground truth motif covered
%                   3: the ground truth occurrence covered
%                   4: the fraction of the ground truth occurrence covered
%                   5: the fraction of the found occurrence covering the
%                      ground truth occurrence
%                   6: beginning of overlapping region between ground truth
%                      and found occurrences
%                   7: end of overlapping region between ground truth
%                      and found occurrences
%
% coveredby         a cell array with an element for each motif in the
%                   ground truth motifs. this element is an array of 7 columns
%                   giving the following information in order:
%                   1: the occurrence of the ground truth covered
%                   2: the found motif covering
%                   3: the found occurrence covering
%                   4: the fraction of the ground truth occurrence covered
%                   5: the fraction of the found occurrence covering the
%                      ground truth occurrence
%                   6: beginning of overlapping region between ground truth
%                      and found occurrences
%                   7: end of overlapping region between ground truth
%                      and found occurrences
%
%
% Notes:
% ======
% * (1-coveringNone-coveringMulti-coveringPartial) gives the fraction of 
%    found motifs that
%    are clearly covering a single motif in the ground truth (these are our
%    best friends!!)
% *  in case the number of learned motifs is zero (locs is empty), all
%    covering* are returned equal to zero except coveringNone which will be
%    equated to 1
%
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, G-SteX: Greedy Stem Extension for
% Free-Length Constrained Motif Discovery, IEA/AIE 2012
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%
if nargin<3 || isempty(minCover2Consider)
    minCover2Consider=0;
end
if nargin<4 || isempty(maxOverlap2Neglect)
    maxOverlap2Neglect=minCover2Consider; %max(0.1,2*minCover2Consider);
end
nTrue=length(locsT);
nLearned=length(locs);

ignoreMultiAndPartialInCoverage=0;

coveringNone=0;
coveringMulti=0;
coveringPartial=0;
badMotif=zeros(nLearned,1);
covered=zeros(nTrue,1);
totals=zeros(nTrue,1);
extras=zeros(nTrue,1);
covering=cell(nLearned,1);
coveredby=cell(nTrue,1);
if nLearned==0
    coveringNone=1;
    return;
end
for i=1:nTrue    
    nt=size(locsT{i},1);
    for kt=1:nt        
        for j=1:nLearned
            no=size(locs{j},1);
            for l=1:no
                [fraction,reverse,covRange]=...
                    coveringFraction(locsT{i}(kt,:),locs{j}(l,:));
                if fraction>minCover2Consider
                    covering{j}=[covering{j};l,i,kt,fraction,reverse,covRange];
                    coveredby{i}=[coveredby{i};kt,j,l,fraction,reverse,covRange];
                end
            end
        end
    end    
end
for i=1:nLearned
    if isempty(covering{i})
        coveringNone=coveringNone+1;
        badMotif(i)=1;
        continue
    end
    if size(covering{i},1)<size(locs{i},1)
        coveringPartial=coveringPartial+1;
        badMotif(i)=2;
        continue
    end
    myMotifs=unique(covering{i}(:,2));    
    if numel(myMotifs)>1
        toremove=zeros(size(myMotifs));
        for j=numel(myMotifs(:)')
            toremove(j)=max(covering{i}(covering{i}(:,2)==myMotifs(j),5))...
                <=maxOverlap2Neglect;
        end        
        myMotifs(toremove==1)=[];
    end
%     if isempty(myMotifs)      % this is commented out because this case
%         coveringNone=coveringNone+1;  % will be found in extras
%         badMotif(i)=1;
%     else
    if numel(myMotifs)>1
        coveringMulti=coveringMulti+1;
        badMotif(i)=3;    
    end
end
coveringMulti=coveringMulti/nLearned;
coveringNone=coveringNone/nLearned;
coveringPartial=coveringPartial/nLearned;

if nargout<3
    return;
end

if ignoreMultiAndPartialInCoverage
    for i=1:nTrue
        nC=size(coveredby{i},1);
        toremove=zeros(nC,1);
        for j=1:nC
            if badMotif(coveredby{i}(j,2))~=0
                toremove(j)=1;
            end
        end
        coveredby{i}(toremove==1,:)=[];
    end
end
for i=1:nTrue
    nt=size(locsT{i},1);        
    for k=1:nt        
        me=zeros(locsT{i}(k,2)-locsT{i}(k,1)+1,1);
        if ~isempty(coveredby{i})
            rangs=coveredby{i}(coveredby{i}(:,1)==k,[6,7]);
            rangs=rangs-locsT{i}(k,1)+1;
            rangs(rangs<1)=0;
            for r=1:size(rangs,1)
                me(rangs(r,1):rangs(r,2))=ones(rangs(r,2)-rangs(r,1)+1,1);
            end
            covered(i)=covered(i)+sum(me);            
        end        
        totals(i)=totals(i)+numel(me);
        if nargout<4
            continue;
        end          
        if ~isempty(coveredby{i})
            mycover=coveredby{i}(coveredby{i}(:,1)==k,[2,3,5]);        
            for r=1:size(mycover,1)
                extras(i)=extras(i)+round((1-mycover(r,3))*  ...
                    (locs{mycover(r,1)}(mycover(r,2),2)- ...
                     locs{mycover(r,1)}(mycover(r,2),1)+1));
            end        
        end
    end    
end
covered=covered./totals;
extras=extras./totals;

end


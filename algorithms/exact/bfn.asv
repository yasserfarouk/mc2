function [locs,...
    dists,executionTime,means,stats]=bfn(x,lrange,lStep,normalization,innerOverlap ...
                        ,nMotifs,outerOverlap,fast)
% applies the MK algorithm using mk++.exe which is a slightly modified version
% of mk_l.exe supplied by Mueen and Keogh.
% parameters are the same as parameters for MK
%
% Inputs:
% =======
% x         The input time series
% lrange    The minimum and maximum expected lengths of motifs. If a
%           single number is passed it is treated as the maximum 
% lStep     The number of points to add from the lowest to highest motif lengths.
%           if not given then it is assumed min(10,(maxLength-minLength)/10)
% normalization  [Optional] The normalization method used:
%            0  NO normalization
%            1  subtract the mean
%            2  subtract the mean and divide by the range 
%            3  zscore (subtract the mean and divide by the standard
%               deviation)            (default)
% innerOverlap = [Optional] Maximum overlap allowed between occurrences in the same
%                motif. Window of Exclusion is calculated as: 
%                (1-innerOverlap)*LENGTH. The algorithm does not calculate 
%                distances between any pairs that start less than W points 
%                apart (default 0.0 allowing no overlap).
% R = [Optional] number of reference points (default 8)
% nMotifs =  [Optional] number of motifs to discover. (default 10)
% outerOverlap = [Optional] fraction of overlap allowed between different
%               motifs (default 0.0=no overlap)
% fast      = [Optional] if nonzero, then best-so-far will be used to prune
%               distance caclulations (default 1)
% Outputs:
% ========
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
% maxSmallDistance  The Euclidean distance between the two occurrences of 
%                   the last motif found
% success   nonzero if the function succeeded. if it fails (due to a
%           problem with mk++.exe) it returns zero


    if ~isvector(x)
        if size(x,1)>size(x,2); x=x'; end;
        x=tspca(x);
    end
    if length(lrange)<2
        lrange=[max(3,floor(0.9*lrange)),lrange];
    end
    if nargin<3
        lStep=max(10,ceil((lrange(2)-lrange(1))/10));   
    end
    if nargin<4
        normalization=3;
    end
    if nargin<5
        innerOverlap=0;
    end    
    if nargin<6
        nMotifs=10;
    end
    if nargin<7
        outerOverlap=0.0;
    end          
    if nargin<8
        fast=1;
    end   
    nMotifLengths=floor(1+(lrange(2)-lrange(1))/lStep);
    locs=[]; means=[]; stats=[];
    success=0;
    x=x(:);    
    current=cd();
    p= mfilename('fullpath');
    p(end-3:end)=[];
    p=sprintf('%s/bin',p);
    cd(p);
    csvwrite(strcat(p,'/tmpMKData.csv'),x);
    str=sprintf('bfn.exe tmpMKData.csv %d %d %d %f %d %d %f %d %d', ...
        length(x),lrange(1),lrange(2),lStep,innerOverlap,nMotifs,outerOverlap,normalization,fast);
    [s,result]=system(str);    
    cd(current);
    if s
        executionTime=inf;
        return;
    end
    success=1;            
    if(nargout>3)
        [locs,dists,executionTime,means,stats]=parseMKoutput(result,x);
    else
        [locs,dists,executionTime]=parseMKoutput(result,x);
    end            
end

function s=square(x)
s=x.*x;
end
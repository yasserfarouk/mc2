function [locs,means,stats,maxSmallDistance,executionTime,success]=bfmd(x,lengths,innerOverlap ...
                        ,nMotifs,outerOverlap,doNormalize, Fast)
% applies the brute force algorithm
%
% Inputs:
% =======
% lengths = if scalar --> motif length
%     if a vector --> motif discovery is applied at all lengths in this
%     vector (e.g. [2:2:10] will apply bruteForce to lengths 2,4,6,8,10)
% innerOverlap = Maximum overlap allowed between occurrences in the same
%                motif. Window of Exclusion is calculated as: 
%                (1-innerOverlap)*LENGTH. The algorithm does not calculate 
%                distances between any pairs that start less than W points 
%                apart (default 0.0 allowing no overlap).
% nMotifs =  [Optional] number of motifs to discover. (default 1)
% outerOverlap = [Optional] fraction of overlap allowed between different
%               motifs (default 0.0=no overlap)
% doNormalize = [Optional] if nonzero the data is normalized by subtacting mean and
%               dividing by std deviation (default 1)
% F =        [Optional] if nonzero then only the maximum not max K motifs are kept during
%           first stage of references (default 1)

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
        x=tspca(x);
    end
    if nargin<3
        innerOverlap=l;
    end    
    if nargin<4
        nMotifs=1;
    end
    if nargin<5
        outerOverlap=0.0;
    end
    if nargin<6
        doNormalize=1;
    end
    if nargin<7
        Fast=1;
    end
    executionTime=0;
    maxSmallDistance=0;
    locs=[]; means=[]; stats=[];
    success=0;
    x=x(:);    
    current=cd();
    p= mfilename('fullpath');
    p(end-4:end)=[];
    cd(p);
    csvwrite(strcat(p,'tmpMKData.csv'),x);
    locs=cell(nMotifs*length(lengths),1);
    nxtLoc=1;
    for l=lengths
        str=sprintf('bruteForce.exe tmpMKData.csv %d %d %f %d %f %d %d', ...
            length(x),l,innerOverlap,nMotifs,outerOverlap,doNormalize, Fast);
        [s,result]=system(str);    
        cd(current);
        if s
            return;
        end
        success=1;
        nums=textscan(result,'%*s %*s %*s %*s %*s %*s %d%*s %d %*s %*s %*s %*s %*s %*s %f');
        eTime=textscan(result(strfind(result,'Time:')+6:end),'%f %*s');
        executionTime=executionTime+eTime{1};
        nums{1}(end,:)=[];      % remove the extra line containing execution time
        nums{2}(end,:)=[];
        nums{3}(end,:)=[];
        nMotifs=length(nums{1});
        if nMotifs<1
            return;
        end
        
        for i=1:nMotifs
            % notice that C arrays are zero based while MATLAB's are 1-based
            locs{nxtLoc}=double([nums{1}(i)+1,nums{1}(i)+l;nums{2}(i)+1,nums{2}(i)+l]);
            [stats,means]=updateStatsAndMean(locs{nxtLoc},x,@square,stats,means,nxtLoc);
            nxtLoc=nxtLoc+1;
        end
        maxSmallDistance=max(maxSmallDistance,max(nums{3}));
    end
end

function s=square(x)
s=x.*x;
end
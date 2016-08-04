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

function [locs,...
    dists,executionTime,means,stats]=mkp(x,l,innerOverlap,R ...
                        ,nMotifs,outerOverlap,doNormalize, Fast)
% applies the MK algorithm using mk++.exe which is a slightly modified version
% of mk_l.exe supplied by Mueen and Keogh.
% parameters are the same as parameters for MK
%
% Inputs:
% =======
% l = motif length
% innerOverlap = Maximum overlap allowed between occurrences in the same
%                motif. Window of Exclusion is calculated as: 
%                (1-innerOverlap)*LENGTH. The algorithm does not calculate 
%                distances between any pairs that start less than W points 
%                apart (default 0.0 allowing no overlap).
% R = [Optional] number of reference points (default 8)
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
        innerOverlap=0;
    end
    if nargin<4
        R=8;
    end
    if nargin<5
        nMotifs=1;
    end
    if nargin<6
        outerOverlap=0.0;
    end
    if nargin<7
        doNormalize=1;
    end
    if nargin<8
        Fast=1;
    end
    locs=[]; means=[]; stats=[];
    success=0;
    x=x(:);    
    current=cd();
    p= mfilename('fullpath');
    p(end-2:end)=[];
    p=sprintf('%s/bin',p);
    cd(p);
    csvwrite(strcat(p,'/tmpMKData.csv'),x);
    str=sprintf('mk+.exe tmpMKData.csv %d %d %f %d %d %f %d %d', ...
        length(x),l,innerOverlap,R,nMotifs,outerOverlap,doNormalize, Fast);
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
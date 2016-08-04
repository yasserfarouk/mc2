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

function [x, motifs,locs, occurs,y]=generatePatterns(nMotifs,occurences,pLen,pRandBefore,randRangeBefore,randLenMax,noiseAmplitude,gNoiseSigma,L,randLocs)
% A function to generate time series with recurring motifs in them
% interleaved with random data
% 
% nMotifs             number of motif patterns used. The first three
%                     patterns are always a sinusoidal, linear, dc (of
%                     course if you request at least 3 patterns)
% occurences          array of nMotifs elements specifying how many
%                     occursrencs of each motif to use
% pLen=[100 100];     minimum and maximum length of each pattern
% pRandBefore=0.3;    probability of having a random part before each pattern
% randRangeBefore=[0 0];    range of random data added before each pattern if any
% randLenMax=2*pLen;  maximum length of random data added before each pattern 
% noiseAmplitude=0;   random purterpation added to each pattern (the actual
%                     noise will be between -noiseAmplitude/2 and noiseA../2
% gNoiseSigma=1;      the variance of a gaussian noise added to the
%                     patterns and the whole signal so the noise will be from N(0,gNoiseSigma)          
% L=0;                minimum length of the timeseries. If zero then it is
%                     not controlled otherwise the timeseries is extended
%                     to this length
% randLocs=0         if 0 then the locations of motifs are distributed
%                    uniformly on the length of the timeseries
%                    if 1 then the locations of motifs are distributed
%                    randomly on the length of the timeseries
%
%output
%======
% x                The time series
% locs             a cell of nMotifs elements specyfing the occursrance
%                  locations of each motif
% occurs            A cell array giving the actual occursrances
%                  occursances{i,j}  =  the occursrance j of motif i
% motifs           a cell array with the tru motifs
% y                a time series representing x without any noise
%
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%

if~exist('nMotifs','var')
    nMotifs=4;
end
if(nMotifs<1)
    error('Cannot make a sequence with no motifs');
end
if~exist('occurences','var')
    occurences=10*ones(1,nMotifs);
end
if length(occurences)<nMotifs
    error('number of elements in occurences must equal nMotifs');
end
if length(occurences)>nMotifs
    warning('only the first nMotifs of occurences will be used');
end
if~exist('pLen','var')
    pLen=[100,100];
end
if~exist('pRandBefore','var')
    pRandBefore=0.3;
end
if~exist('randRangeBefore','var')
    randRangeBefore=[-0.5 0.5];
end

if~exist('randLenMax','var')
    randLenMax=2*pLen;
end
if~exist('noiseAmplitude','var')
    noiseAmplitude=0;
end
if~exist('gNoiseSigma','var')
    gNoiseSigma=0;
end

if~exist('L','var')
    L=0;
end

if~exist('randLocs','var')
    randLocs=0;
end
smoothingFactor=2;
        
motifs=cell(nMotifs,1);
lengths=zeros(nMotifs,1);
done =0;
if(rand(1,1)>0.5)
    peakpeak=randRangeBefore(2)-randRangeBefore(1);
    if peakpeak<0.0001
        peakpeak=1;
    end
    if(nMotifs>0)
        l=round((pLen(2)-pLen(1))*rand(1,1)+pLen(1));
        motifs{1,1}=.5*peakpeak*sin(0:2*pi/l:2*pi*(1-1/l)); % first pattern is a sin
        done=done+1;
    end
end
% if(nMotifs>1)
%     l=round((pLen(2)-pLen(1))*rand(1,1)+pLen(1));
%     motifs{2,1}=1*ones(1,l);            % second pattern is linear
%     done=done+1;
% end
% if(nMotifs>2)
%     l=round((pLen(2)-pLen(1))*rand(1,1)+pLen(1));
%     motifs{3,1}=-1:2/l:(1-1/l);            % third pattern is linear
%     done=done+1;
% end
if nMotifs>done    
    for i=done+1:nMotifs
        l=round((pLen(2)-pLen(1))*rand(1,1)+pLen(1));
        smoothingFactor=max(2,min(floor(l/6),max(7,floor(l/12))));
        motifs{i,1}=smooth(rand(1,l)-0.5,smoothingFactor);             % pattern is randomly selected
        motifs{i,1}=(.5/max(motifs{i,1})).*motifs{i,1};
        motifs{i,1}=motifs{i,1}(:)';
    end
end
for i=1:nMotifs
    lengths(i)=numel(motifs{i,1});
end
locs=cell(nMotifs,1);
%lengths=cell(nMotifs,1);
occurs=cell(nMotifs,1);
for i=1:nMotifs
    locs{i,1}=zeros(occurences(i),2);
  %  lengths{i,1}=zeros(1,occurences(i));
    occurs{i,1}=cell(occurences(i),1);
end
x=[];
y=[];
rem=occurences;
%will be locations of all occurences of all motifs
os=zeros(sum(occurences),1);  
osDone=0;
for i=1:nMotifs
    os(osDone+1:osDone+occurences(i))=i.*ones(occurences(i),1);
    osDone=osDone+occurences(i);
end
nTotal=numel(os);
os=os(randperm(nTotal));

for osDone=1:nTotal;
    if(rand(1,1)<=pRandBefore)
        rl=round(randLenMax*rand(1,1));
        if(rl>0)
            randSeg=(randRangeBefore(2)-randRangeBefore(1))*rand(1,rl)+randRangeBefore(1);
            randSeg=smooth(randSeg,smoothingFactor)';
            %randSeg=(.5/max(randSeg)).*randSeg;
            x=[x randSeg];
            y=[y zeros(1,rl)];
        end
    end    
    p=os(osDone);     
    pattern=motifs{p,1};    
    pp=pattern+noiseAmplitude.*(rand(1,length(pattern))-0.5*ones(1,length(pattern)));
    current=occurences(p)-rem(p)+1;
    locs{p,1}(current,1)=length(x)+1;
    locs{p,1}(current,2)=length(x)+length(pp);
    %lengths{p,1}(current)=length(pp);    
    occurs{p,1}{current}=pp;
    x=[x pp];
    y=[y pattern];
    rem(p)=rem(p)-1;    
end

if(rand(1,1)<=pRandBefore)
    rl=round(randLenMax*rand(1,1));
    if(rl>0)
        randSeg=(randRangeBefore(2)-randRangeBefore(1))*rand(1,rl)+randRangeBefore(1);
        randSeg=smooth(randSeg,smoothingFactor)';
        x=[x randSeg];
        y=[y zeros(1,rl)];
    end
end


if (L~=0)
    nAll=numel(x);
    if nAll>L
        warning('MATLAB:length','You specify %d length but we needed %d points',L,nAll);
    elseif nAll<L
       before=floor((L-nAll).*rand(1,1));
       after=L-nAll-before;
       if before>0
            z=(randRangeBefore(2)-randRangeBefore(1))*rand(1,before)+randRangeBefore(1); 
            z=smooth(z,smoothingFactor)';
            y=[zeros(size(z)) y];
            x=[z x];
            for i=1:nMotifs
                nz=numel(z);
                locs{i}=locs{i}+nz*ones(occurences(i),2);                
            end
       end
       if after>0
           z=(randRangeBefore(2)-randRangeBefore(1))*rand(1,after)+randRangeBefore(1); 
           z=smooth(z,smoothingFactor)';
           y=[y zeros(size(z))];
           x=[x z];
       end
    end
end

x=x+gNoiseSigma*randn(size(x));
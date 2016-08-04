function [stats,delayStats]=cpquality(foundLocs,trueLocs,T,nPoints,useBalancedStats,localizationFun,varargin)
% finds the quality of change points in fouldLocs using trueLocs as true change point
% locaitons. 
%
% foundLocs if foundLocs is a cell of m*1 elements then each element is an array
%           representing locations of change according to one algorithm
%           if foundLocs is n*m array then it represents m outputs of m algorithms with a time series of
%           length n. This means that each algorithm's output must be in a
%           column not a row. The localization function and threshold will
%           be used to create a cell array of locations to be used in the
%           evaluation
% trueLocs  true change points in a vector or cell array
% T         the length of the time series for which these locations where
%           found. This is important to find the number of negatives
% nPoints   width around change points that is assumed OK. If a negative
%           number it is selected as distance to next true chagne point
%           divided by -nPoints
% localizationFun   The function used to localize changes (default is
%                   findLocsThMax. Pass -1 if you want to use the default
% varargin The parameters to be passed to the localization function
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


% check input
if(nargin<4)
    nPoints=0;
end
if(nargin<5)
    useBalancedStats=1;
end
if(nargin<6 || ~ishandle(localizationFun))
    localizationFun=@findLocsTh;
end

if(~iscell(foundLocs))
    if(size(foundLocs,1)==1) && (size(foundLocs,2)~=1)
        foundLocs=foundLocs';
    end
    x=foundLocs;
    m=size(x,2);
    foundLocs=cell(m,1);
    for i=1:m
        foundLocs{i}=localizationFun(x(:,i),deal(varargin{:}));
    end
end

if(iscell(trueLocs))
    trueLocs=cell2mat(trueLocs);
end

if ~isvector(trueLocs)
    error('trueLocs must be a vector');
end

m=numel(foundLocs);

% initialize
stats=[];
delayStats=[];

for c=1:m       % for each column of foundLocs        
    [tp,fp,tn,fn,dMean,dStd,dMax,dMin]=mychangeQuality(foundLocs{c},trueLocs,nPoints,T,useBalancedStats);
    
    if((tp+fn)<eps) || (tp+fn)<eps || (tn+fp)< eps || (tn+fn)<eps
        mcc=(tp*tn-fp*fn);
    else
        mcc=(tp*tn-fp*fn)/(sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)));
    end

    p=fn+tp;
    n=(fp + tn);
    if p==0;
        recall=0;
    else
        recall=tp / p;
    end
    if (tp + fp)==0
        precision=0;
    else
        precision= tp / (tp + fp);
    end
    if (fp + tp)==0
        fdr=1;
    else
        fdr=fp/(fp + tp);
    end
    if  (tn + fn)==0
        npv=0;
    else
        npv=tn / (tn + fn);
    end
    spec=tn / n;
    if (precision+recall)==0
        F=0;
    else
        F=2*precision*recall/(precision+recall);
    end
    stats=[stats;struct('fp',fp,'fn',fn,'tp',tp,'tn',tn...
                ,'mcc',mcc...
                ,'f1',F,...
                'precision',precision,'recall',recall,...
                'accuracy', (tp + tn) / (p + n),'specificity',spec, ...
                 'sensitivity',recall ...
                ,'fpr', fp / n,'fnr',fn/p ...
                ,'ppv',precision,'npv',npv...
                ,'fdr', fdr)];
    delayStats=[delayStats;struct('mean',dMean,'StdDev',dStd,'max',dMax,'min',dMin)];
end


end

% foundLocs here MUST be a vector and assumes that its maximum is 1
function [tp,fp,tn,fn,dMean,dStd,dMax,dMin]=mychangeQuality(myLocs,trueLocs,nPoints,T,useBalancedStats)
tp=0; 
trueLocs=sort(trueLocs); trueLocs=trueLocs(:); n=numel(trueLocs);
myLocs=sort(myLocs); myLocs=myLocs(:);
ntrue=numel(trueLocs); nmy=numel(myLocs);
used=boolean(zeros(nmy,1));
found=boolean(zeros(ntrue,1));
delays=[];
if(ntrue==0)
    tp=0;
elseif(nmy==0)
    tp=0;
else
    % find the location of nearest change to every change point
    [delays,locs]=getDelays(trueLocs,myLocs);

    %find the number of true positives
    if(nPoints>=0)
        %tp=numel(delays(abs(delays)<=nPoints));    
        for i=1:n        
            if (abs(delays(i))<=nPoints)     
                if(~used(locs(i)))
                    tp=tp+1; used(locs(i))=1; found(i)=1;
                end
            end
        end
    else    
        difToBefore=(trueLocs(1:n)-[1;trueLocs(1:n-1)])./-nPoints;
        difToNext=([trueLocs(2:n);T]-trueLocs(1:n))./-nPoints;
        for i=1:n        
            if (delays(i)<0 && (-delays(i))<difToBefore(i)) || ...
               (delays(i)>=0 && (delays(i))<difToNext(i))     
                if(~used(locs(i)))
                    tp=tp+1; used(locs(i))=1; found(i)=1;
                end
            end
        end    
    end
end
fn=ntrue-tp;
fp=nmy-sum(used);
if(useBalancedStats)    % make the number of negatives roughly equal to positives
    negatives=max(fn,round((T-nmy)*(double(ntrue)/double(T))));
    tn=negatives-fn;
else
    negatives=T-nmy;
    tn=negatives-fn;
end

if isempty(delays)
    dMean=inf;
    dStd=inf;
    dMax=inf;
    dMin=-inf;
else
    if isempty(delays(found==1))
        dMean=inf;
        dStd=inf;
        dMax=inf;
        dMin=-inf;
    else
        delays=double(delays);
        dMean=mean(delays(found==1)); if isnan(dMean); dMean=inf; end;
        dMax=max(delays(found==1)); if isempty(dMax); dMax=inf; end;
        dMin=min(delays(found==1)); if isempty(dMin); dMin=-inf; end;  
        if(dMax==dMin)
            dStd=0;
        else
            dStd=std(delays(found==1));   if isnan(dStd); dStd=inf; end;              
        end
    end
end

end

function [ds,locs]=getDelays(x,y)%,rmBottom,rmTop)
% calculates the delays from changes in x to changes in y
    ds=[];    
    locs=zeros(size(x));
    n=numel(x); m=numel(y);
    if n==0 && m==0 
        return;
    end    
    for i=1:n
        fnd=0;
        for j=1:m
            if(y(j)>=x(i))
                d=y(j)-x(i);
                if(j>1)
                    db=x(i)-y(j-1);  locs(i)=j-1;
                else
                    db=inf;
                end
                if(db<d)
                    ds=[ds;-db];  locs(i)=j-1;
                else
                    ds=[ds;d];  locs(i)=j;
                end
                fnd=fnd+1;
                break;
            end            
        end
        if(~fnd)
            ds=[ds;y(end)-x(i)];
            locs(i)=length(y);
        end
    end
end
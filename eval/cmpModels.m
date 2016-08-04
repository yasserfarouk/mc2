function [stats,timeStats] = cmpModels(t,m,calcSelfLoops)
% finds error statistics for the model m compared with the model t (both
% are causal models (see produceMulti for their meaning)
nt=size(t,1);
nm=size(m,1);
if(nt~=nm)
    error('the two models must contain the same number of processes/timeseries');
end
if nargin<3
    calcSelfLoops=0;
end
tp=0; fp=0; tn=0; fn=0;
meanDiffs=[];
stdDiffs=[];
tMeans=[];
tStd=[];
pChgs=[];
for i=1:nt
    if isempty(t{i,1})
        % there can be no true positive all positive here are false
        if isempty(m{i,1})
            % m gets all negatives in this row right
            if(calcSelfLoops);   tn=tn+nt; else  tn=tn+nt-1;  end;            
        else
            % m gets some false positives
            fp=fp+size(m{i,1},1);
            if(calcSelfLoops);   tn=tn+nt-size(m{i,1},1); else  tn=tn+nt-1-size(m{i,1},1);  end;
        end
    else
        % there are some positives 
        if isempty(m{i,1})
            % m lost all positives in this row
            fn=fn+size(t{i,1},1);
            if(calcSelfLoops);   tn=tn+nt-size(t{i,1},1); else  tn=tn+nt-1-size(t{i,1},1);  end;                        
        else
            % both m ant t got something . we will check them one by one
            for j=1:nt
                if (~calcSelfLoops) && j==i
                    continue;
                end
                tind=find(t{i,1}(:,1)==j,1);
                if(isempty(tind))
                    % j is not a positive so we should not find it in m
                    mind=find(m{i,1}(:,1)==j,1);
                    if(isempty(mind))
                        % we did not find the negative correct
                        tn=tn+1;
                    else
                        fp=fp+1;
                    end
                else
                    % this is a positive so we should find it in m
                    mind=find(m{i,1}(:,1)==j,1);
                    if(isempty(mind))
                        % we did not find this positive
                        fn=fn+1;
                    else
                        tp=tp+1;
                        % here we update timeStats
                        if t{i,1}(2)>-eps && m{i,1}(2)>-eps
                            meanDiffs=[meanDiffs;(m{i,1}(2)-t{i,1}(2))];
                            stdDiffs=[stdDiffs;(m{i,1}(3)-t{i,1}(3))];
                            tMeans=[tMeans;t{i,1}(2)];
                            tStd=[tStd;t{i,1}(3)];
                            pChgs=[pChgs;t{i,1}(4)];
                        end
                    end
                end
            end            
        end
    end
end
if calcSelfLoops
    assert(nm*nm==(fp+fn+tp+tn),'not everyting is accounted for');
else
    assert((nm*nm-nm)==(fp+fn+tp+tn),'not everyting is accounted for');
end


if((tp+fn)<eps) || (tp+fn)<eps || (tn+fp)< eps || (tn+fn)<eps
    mcc=(tp*tn-fp*fn);
else
    mcc=(tp*tn-fp*fn)/(sqrt((tp+fn)*(tp+fn)*(tn+fp)*(tn+fn)));
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
stats=struct('fp',fp,'fn',fn,'tp',tp,'tn',tn...
            ,'mcc',mcc...
            ,'f1',F,...
            'precision',precision,'recall',recall,...
            'fpr', fp / n,'fnr',fn/p,'accuracy', (tp + tn) / (p + n),...
            'specificity',spec,'ppv',precision,'npv',npv...
            ,'fdr', fdr,'sensitivity',recall);
timeStats=struct('trueMeans',tMeans,'trueStdDev',tStd,...
                    'meanDiffs',meanDiffs,'stdDiffs',stdDiffs,....
                    'changeProbability',pChgs);
end


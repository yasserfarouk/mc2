function stats = mdq( locs,locsT,T)
%Finds the confusion matrix statistics of the motif discovered. Results are
%aggregated over all motifs
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

nTrue=length(locsT);
nLearned=length(locs);

truex=zeros(T,1);
learnedx=zeros(T,1);                
for j=1:nLearned
    no=size(locs{j},1);
    for l=1:no
        if locs{j}(l,2)<=T && locs{j}(l,1)>0
            learnedx(locs{j}(l,1):locs{j}(l,2))=1;
        end
    end
end
for i=1:nTrue
    nt=size(locsT{i},1);
    for kt=1:nt
        truex(locsT{i}(kt,1):locsT{i}(kt,2))=1;
    end
end
tp=sum(truex & learnedx)./T;
fp=sum(~truex & learnedx)./T;
tn=sum(~truex & ~learnedx)./T;
fn=sum(truex & ~learnedx)./T;

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
stats=struct('fp',fp,'fn',fn,'tp',tp,'tn',tn...
        ,'mcc',mcc...
        ,'f1',F,...
        'precision',precision,'recall',recall,...
        'accuracy', (tp + tn) / (p + n),'specificity',spec, ...
         'sensitivity',recall ...
        ,'fpr', fp / n,'fnr',fn/p ...
        ,'ppv',precision,'npv',npv...
        ,'fdr', fdr);
end


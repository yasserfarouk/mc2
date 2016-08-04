function statsM = mdqM( locs,locsT,T,minOverlapFraction)
%Finds the confusion matrix statistics of the motif discovered. Results are
% calculated for each motif alone
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
fp=zeros(nTrue,1);
tp=zeros(nTrue,1);
fn=zeros(nTrue,1);
tn=zeros(nTrue,1);
for i=1:nTrue
    truex=zeros(T,1);
    learnedx=zeros(T,1);
    nt=size(locsT{i},1);
    for kt=1:nt
        truex(locsT{i}(kt,1):locsT{i}(kt,2))=1;
        for j=1:nLearned
            no=size(locs{j},1);
            for l=1:no
                [isover,fraction]=testOverlapping(locsT{i}(kt,:),...
                    locs{j}(l,:));
                if (minOverlapFraction==0 && isover) || fraction>=minOverlapFraction
                    if locs{j}(l,2)<=T && locs{j}(l,1)>0
                        learnedx(locs{j}(l,1):locs{j}(l,2))=1;
                    end
                end
            end
        end
    end
    tp(i)=sum(truex & learnedx)/T;
    fp(i)=sum(~truex & learnedx)/T;
    tn(i)=sum(~truex & ~learnedx)/T;
    fn(i)=sum(truex & ~learnedx)/T;    
end

p=fn+tp;
n=(fp + tn);
spec=tn ./ n;

mcc=zeros(size(tp));
recall=zeros(size(tp));
precision=zeros(size(tp));
fdr=zeros(size(tp));
npv=zeros(size(tp));
F=zeros(size(tp));
for i=1:nTrue
    if((tp(i)+fn(i))<eps) || (tp(i)+fn(i))<eps || (tn(i)+fp(i))< eps || (tn(i)+fn(i))<eps
        mcc(i)=(tp(i)*tn(i)-fp(i)*fn(i));
    else
        mcc(i)=(tp(i)*tn(i)-fp(i)*fn(i))/(sqrt((tp(i)+fp(i))*(tp(i)+fn(i))*(tn(i)+fp(i))*(tn(i)+fn(i))));
    end
    if p(i)==0;
        recall(i)=0;
    else
        recall(i)=tp(i) / p(i);
    end
    if (tp(i) + fp(i))==0
        precision(i)=0;
    else
        precision(i)= tp(i) / (tp(i) + fp(i));
    end
    if (fp(i) + tp(i))==0
        fdr(i)=1;
    else
        fdr(i)=fp(i)/(fp(i) + tp(i));
    end
    if  (tn(i) + fn(i))==0
        npv(i)=0;
    else
        npv(i)=tn(i) / (tn(i) + fn(i));
    end
    if (precision(i)+recall(i))==0
        F(i)=0;
    else
        F(i)=(2*precision(i)*recall(i))/(precision(i)+recall(i));
    end
end

statsM=struct('fp',fp,'fn',fn,'tp',tp,'tn',tn...
        ,'mcc',mcc...
        ,'f1',F,...
        'precision',precision,'recall',recall,...
        'accuracy', (tp + tn) ./ (p + n),'specificity',spec, ...
         'sensitivity',recall ...
        ,'fpr', fp ./ n,'fnr',fn./p ...
        ,'ppv',precision,'npv',npv...
        ,'fdr', fdr);
end


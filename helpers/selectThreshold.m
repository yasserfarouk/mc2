function th=selectThreshold(x,L,distFun)
if nargin<3 || isempty(distFun)
    distFun=@dE;
end
T=numel(x);

nSamples=T;
fraction=1/30.0;
minSeparation=L;
nTrials=3;
tmpLoc=randi(T-L+1,nSamples,1);
tmpLoc2=randi(T-L+1,nSamples,1);
indx=find(abs(tmpLoc-tmpLoc2)<minSeparation);
for i=1:nTrials
    
    if isempty(indx)
        break;
    end
    nBadSamples=numel(indx);
    tmpLoc(indx)=randi(T-L+1,nBadSamples,1);
    tmpLoc2(indx)=randi(T-L+1,nBadSamples,1);
    indx=find(abs(tmpLoc-tmpLoc2)<minSeparation);
end

tmpLoc(indx)=[];
tmpLoc2(indx)=[];
nSamples=numel(tmpLoc);
d=zeros(nSamples);
for i=1:nSamples
    d(i)=dE(x(tmpLoc(i):tmpLoc(i)+L-1),x(tmpLoc2(i):tmpLoc2(i)+L-1));
end
d=sort(d);
th=d(ceil(fraction*T));
end
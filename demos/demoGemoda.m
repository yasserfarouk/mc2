x=[1,2,3,4,5,6,1,2,3,5,4,6]';
locs=gemoda(x,3,2,2,@hamming,true);
showMDResults(x,locs);

locs=gemoda(x,3,1,2,@hamming,true);
showMDResults(x,locs);



N=200;
L=50;
T=1000;
range=[-5,5];    
scale2range=1;
noiseLevel=0.00;
gNoiseSigmaInitial=0.0;       % initial noise
pOutlier=0.0;
minChangeDistance=50;
pchange=0.01;

ps=[5,10,30];
qs=[0];

[xbase,locsTbase,tbase,xTruebase]=produceSingle(T,'pchange',pchange,...
    'minChangeDistance',minChangeDistance,'noiseLevel',noiseLevel,...
    'poutlier',pOutlier,'noiseSigma',gNoiseSigmaInitial,'range',range...
    ,'scale2range',scale2range);

Lmin=L/2;
tmpLoc=randi(T-Lmin+1,T);
tmpLoc2=randi(T-Lmin+1,T);
d=zeros(T);
for i=1:T
    d(i)=dE(xbase(tmpLoc(i):tmpLoc(i)+Lmin-1),xbase(tmpLoc2(i):tmpLoc2(i)+Lmin-1));
end
d=sort(d);
th=d(floor(T/30));

[locs]=gemoda(xbase,40,th,2,[],[],[],[],40);
showMDResults(xbase,locs);
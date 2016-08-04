nMotifsToUse=2;
nOccursMin=6;
nOccursMax=10;
minMotif=30;
maxMotif=60;
gSigma=0.00;

th=0.1;
nPointsPerSymbol=5;
nSymbolsPerWord=3;
alphabetSize=15;
Lmin=30;
Lmax=60;

rng('default');

[x,motifsT,locsT,occursT]=generatePatterns(...
            nMotifsToUse,randi([nOccursMin,nOccursMax],1,nMotifsToUse)...
            ,[minMotif,maxMotif],1,[-0.5,0.5]...
            ,4*maxMotif,0,gSigma);
T=numel(x);
remaining=(1+ceil(T/maxMotif))*maxMotif-T;
x=[x (rand(1,remaining)-0.5)]';

c=cpd(x,floor(Lmin/3),Lmin-floor(Lmin/3));
candLocs=randp(c,floor(T/20));
locs=mcfull(x,[Lmin,Lmax],candLocs);
locs=mdRemoveShortMotifs(locs,Lmin);
showMDResultsWithGT(x,locsT,locs);
 
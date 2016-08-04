% Sets variables to call testcmd to find the effect of noise level
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

doFinalizations=0;

generateNew=1;
generateNewNoise=1;
generateNewLocNoise=1;
generateNewCPDLocs=1;
generateNewDGRResults=1;
applyCMD=1;
applyMDFInalization=doFinalizations;
generateNewStats=1;
generateNewFinalizationStats=1;
trySlidingDuringFinalExtension=0;
showResults=1;
delayList=0:5:25;
nSignals=2;
nMotifsToUse=2;
nOccursMax=5;
nOccursMin=3;
minMotif=60; 
maxMotif=120; 
M=10; 
K=12;
strLen=minMotif;
nSymbols=10;
nAlphabit=20;
nIterations=ceil(.7*nSymbols);
fractionKept=0.7;
noises=[0,0.0001,0.001,0.01,0.1];       
        % this is the noise to be added to the signal itself. 
        
locNoises=[0];%[0,0.000001,0.00001,0.001,0.1];
        %    this is the noise to be added to the localization 
        %   (a fraction of  the motif length). This is only effective if
        %   usersst==0
        %   if you use a negative number rsst will be used for this
        %   locNoise value. For example locNoises=[0,.1,-1] will calculate
        %   using localization noises of 0 and .1 and then using rsst


usersst=1;

algNames={'Projections';'MK++';'G-SteXS';'G-SteXB';'MCFull';'Shift Density'};
algFuns={@projections;@mkppC;@dgr;@dgr;@mcfull;@sdCMD};
algUsesCP=[0         ,0     ,1   ,1   ,1      ,1];
algParams={...
        {[]};
        {[]};
        {'tryextendingstemsslowly',1,'usebisection',0};...
        {'tryextendingstemsslowly',0,'usebisection',1};
        {[]};{[]}...
        };
    
testcmd
if showResults
    displayTestResults;
end
    

% Sets the required variables for testcmd to test the effect of
% localization noise
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
nMotifsToUse=3;
nOccursMax=5;
nOccursMin=3;
minMotif=25; 
maxMotif=100; 
minMotif2Gen=50; 
maxMotif2Gen=75; 
nMotifsBetween=2;
nSignals=5;
% notice that the length of the time series will be at least
% nOccursMax*nMOtifsToUse*maxMotif*(1+nMotifsBetween)
doFinalizations=0;

noises=[0.1];       
        % this is the noise to be added to the signal itself. 
nMAX=nOccursMax*nMotifsToUse*maxMotif*(1+nMotifsBetween);        
lengths=[1,2,4];
        %    this is the noise to be added to the localization 
        %   (a fraction of  the motif length). This is only effective if
        %   usersst==0
        %   if you use a negative number rsst will be used for this
        %   locNoise value. For example locNoises=[0,.1,-1] will calculate
        %   using localization noises of 0 and .1 and then using rsst


generateNew=1;
generateNewNoise=1;
generateNewCPDLocs=1;
generateNewDGRResults=1;
applyCMD=1;
applyMDFInalization=0;
generateNewStats=1;
generateNewFinalizationStats=0;
trySlidingDuringFinalExtension=0;
showResults=1;
delayList=0:5:25;


M=10; 
K=12;
strLen=minMotif;
nSymbols=10;
nAlphabit=20;
nIterations=ceil(.7*nSymbols);
fractionKept=0.7;

algNames={'MK--';'MK++';'sMD';'MOEN'};
algFuns={@mkmmC;@mkppC;@mnC;@moenC};
algUsesCP=[0     ,0   ,0   ,0];
algParams={...        
        {[]};
        {[]};...        
        {[]};
        {[]}...
        };
testcmdL
if showResults
    displayTestResultsLengths;
end
    

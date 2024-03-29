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

doFinalizations=0;

noises=[0.05];       
        % this is the noise to be added to the signal itself. 
        
locNoises=[-1];
        %    this is the noise to be added to the localization 
        %   (a fraction of  the motif length). This is only effective if
        %   usersst==0
        %   if you use a negative number rsst will be used for this
        %   locNoise value. For example locNoises=[0,.1,-1] will calculate
        %   using localization noises of 0 and .1 and then using rsst


usersst=1;
generateNew=1;
generateNewNoise=0;
generateNewLocNoise=0;
generateNewCPDLocs=0;
generateNewDGRResults=1;
applyCMD=0;
applyMDFInalization=0;
generateNewStats=0;
generateNewFinalizationStats=1;
trySlidingDuringFinalExtension=0;
showResults=1;
delayList=0:5:25;
nSignals=2;
nMotifsToUse=2;
nOccursMax=6;
nOccursMin=3;
minMotif=60; 
maxMotif=60; 
M=10; 
K=12;
strLen=minMotif;
nSymbols=10;
nAlphabit=20;
nIterations=ceil(.7*nSymbols);
fractionKept=0.7;

algNames={'GEMODA';'MK++';'Shift Density'};
algFuns={@gemoda;@mkppC;@sdCMD'};
algUsesCP=[0         ,0   ,1];
algParams={...
        {[]};
        {[]};
        {'tryextendingstemsslowly',1,'usebisection',0};...
        {'tryextendingstemsslowly',0,'usebisection',1};
        {[]};{[]}...
        };
testcmd
if showResults
    displayTestResultsLNoises;
end
    

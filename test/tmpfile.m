for s=1:nSignals
for l=1:nLengths                
    [y,locsT{s},~]=embedInside(signals{s,1,1},locsT{s},signals{s,1,1},lengths(l),[-0.5,0.5],0,true,true);
    signals{s,1+1,l}=y;
end
end

save('./Data/signals.mat','signals','motifsT','locsT','occursT','M'...
        ,'K','minMotif','maxMotif','strLen','strings','nSymbols','nIterations','fractionKept');     
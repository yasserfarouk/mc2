function [figs]=showMDResultsWithGT(x,locsT,locs,useMotifRaising)
% displays results of MD with the ground truth (or two MD algorithms
if ~exist('useMotifRaising','var') || isempty(useMotifRaising)
    useMotifRaising=false;
end

figs=figure;
%maximize(figs);
subplot(2,1,1);
showMDResultsOnData(x,locsT,useMotifRaising);
ylabel('Ground Truth');
subplot(2,1,2);
showMDResultsOnData(x,locs,useMotifRaising);
ylabel('Discovered');
end

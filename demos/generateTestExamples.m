%rng(pi);
clear all;
close all;
nModels=20;
T=3000;
nProcesses=[4,8];
noiseSigma=.0;
delayRange=[5,12];
stdDelayRange=[0,0];
signalRange=[-10,10];
scale2range=0;
pchangefree=0.05;
percentAR=.25;
smooth=0;


% 1. Generate data randomly
Ms=cell(nModels,1);
Es=cell(nModels,1);
xs=cell(nModels,1);
ndims=cell(nModels,1);
chgLocs=cell(nModels,1);
chgss=cell(nModels,1);
for i=1:nModels
    fprintf('Model number %d of %d ... ',i,nModels);
    np=ceil(nProcesses(1)+(nProcesses(2)-nProcesses(1))*rand(1,1));
    [Ms{i},Es{i}]=generateCausalGraph(np,'sr',stdDelayRange,'meanrange',delayRange); 
    [xs{i},ndims{i},chgLocs{i},chgss{i}]=produceMulti(Ms{i},T,'range',signalRange,'scale2range',scale2range,'smooth',smooth,'par',percentAR);    
    xs{i}=xs{i}+noiseSigma.*randn(size(xs{i}));
    fprintf('Done\n');
end

mkdir('Data');
save('./Data/ccausalitysynthetic.mat');
fprintf('All DONE\n');
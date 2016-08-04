% Tests different variants of cpd
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, On Comparing SSA-based Change Point
% Discovery Algorithms, IEEE SII 2011, 938-945 
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%
displayResults=1;
displaySeriesNum=5;
showErrorBars=0;
saveFigs=0;
calcESRLoc=0;
calcESRN=0;

LoadSaved=0;
    LoadSavedParams=0;
    LoadSavedSignal=0;
    LoadSavedChanges=0;
    LoadSavedChangeLocs=0;
    LoadSavedMetrics=0;
        LoadSavedTraditionalMetrics=0;
        LoadSavedCPSP=0;
        LoadSavedESR=0;
        LoadSavedESRLoc=0;
        LoadSavedESRN=0;
    

if LoadSaved
    LoadSavedParams=1;
    LoadSavedSignal=1;
    LoadSavedChanges=1;
    LoadSavedChangeLocs=1;
    LoadSavedMetrics=1;        
end
if LoadSavedMetrics
    LoadSavedTraditionalMetrics=1;
    LoadSavedCPSP=1;
    LoadSavedESR=1;
    LoadSavedESRLoc=1;
    LoadSavedESRN=1;
end

if LoadSavedParams    
    load('./Data/params.mat');    
else
    useControlAlgorithms=0;
    nSeriesBase=2;
    useBalancedStats=0;
    T=3000;
    range=[-5,5];    
    scale2range=1;
    noiseLevel=0.00;
    gNoiseSigmaInitial=0.00;       % initial noise
    gNoiseSigmas=[0:0.5:2.5];
    uniformNoises=[0];
    pOutlier=0.0;
    minChangeDistance=50;
    pchange=0.01;
    w=max([3,min(10,floor(minChangeDistance/2)-2)]);
    nwdiff=max([1,ceil(w/2)]);
    n=w+nwdiff;
    points=0:w+n-1;    
    thresholds=[];
    for i=-20:.5:-1
        thresholds=[thresholds,10^i];
    end
    thresholds=[thresholds,.105:.05:1];%[1e-12,1e-6,1e-3,1e-2,1e-1];
    
    findLocsFun=@findLocsTh;
    firstW=-.5;        % see findLocsTh for the meaning of this parameter
    
    statNames={'False Postitives','False Negatives','True Postitives', ...
            'True Negatives', 'MCC', 'F', 'Precision', 'Recall', 'Accuracy', ...
            'Specificity', 'Sensitivity'};
    

    %'lMax'      , [-1,3]; ...
    %'doPost'    , [0,1]; ...    

    %npts4h:     0 ==> l is calculated everytime we calculate a hankel matrix
    %           -1 ==> l is calculated once in the very beginning using the
    %                  whole time series in construction of the Hankel matrix
    %calchankelevely:   0 ==> The hankel matrix for the past is calculated only
    %                         once using the whole time series
    %                   1 ==> It is calculated everytime inside the loop

    rsstParams={ ...
         %{'sst3defaults',1};
         {'sst3defaults',1}; ... % this is MZ algorithm
         {'rsstdefaults',1}; ...   % RSST
         
%          {'sst3defaults',1,'dfunname','distAvgVHDist','normalize',1,'calchankelevery',1,'npts4h',-1}; ...
%          {'sst3defaults',1,'dfunname','distAvgVHDist','normalize',1,'calchankelevery',0,'npts4h',0}; ...
%          {'sst3defaults',1,'dfunname','distAvgVHDist','normalize',1,'calchankelevery',0,'npts4h',0,'symnorm',0}; ...
%          {'sst3defaults',1,'dfunname','distAvgVHDist','normalize',1,'calchankelevery',0,'npts4h',-1}; ...                  
%          {'sst3defaults',1,'dfunname','distAvgVHDist','normalize',0,'calchankelevery',1,'npts4h',0}; ... 
%          {'sst3defaults',1,'dfunname','distAvgVHDist','normalize',0,'calchankelevery',1,'npts4h',-1}; ...
%          {'sst3defaults',1,'dfunname','distAvgVHDist','normalize',0,'calchankelevery',0,'npts4h',0}; ...
%          {'sst3defaults',1,'dfunname','distAvgVHDist','normalize',0,'calchankelevery',0,'npts4h',-1}; ...
         

%          {'sst3defaults',1,'dfunname','distAngBetweenSubspaces','normalize',1,'calchankelevery',0,'npts4h',0}; ...
%          {'sst3defaults',1,'dfunname','distAngBetweenSubspaces','normalize',1,'calchankelevery',0,'npts4h',-1}; ...
%          {'sst3defaults',1,'dfunname','distAngBetweenSubspaces','normalize',1,'calchankelevery',0,'npts4h',0,'symnorm',0}; ...
%          {'sst3defaults',1,'dfunname','distAngBetweenSubspaces','normalize',1,'calchankelevery',1,'npts4h',0}; ...
%          {'sst3defaults',1,'dfunname','distAngBetweenSubspaces','normalize',1,'calchankelevery',1,'npts4h',-1}; ...
%          {'sst3defaults',1,'dfunname','distAngBetweenSubspaces','normalize',0,'calchankelevery',0,'npts4h',0}; ...
%          {'sst3defaults',1,'dfunname','distAngBetweenSubspaces','normalize',0,'calchankelevery',0,'npts4h',-1}; ...
%          {'sst3defaults',1,'dfunname','distAngBetweenSubspaces','normalize',0,'calchankelevery',1,'npts4h',0}; ...
%          {'sst3defaults',1,'dfunname','distAngBetweenSubspaces','normalize',0,'calchankelevery',1,'npts4h',-1}; ...

%          {'sst3defaults',1,'dfunname','distAvgEigDist','normalize',1,'calchankelevery',0,'npts4h',0}; ...
%          {'sst3defaults',1,'dfunname','distAvgEigDist','normalize',1,'calchankelevery',0,'npts4h',-1}; ...
         {'sst3defaults',1,'dfunname','distAvgEigDist','normalize',1,'calchankelevery',0,'npts4h',0,'symnorm',0}; ...
         {'dfunname','distAvgEigDist','normalize',1,'calchankelevery',0,'npts4h',0,'symnorm',0}; ...
%          {'sst3defaults',1,'dfunname','distAvgEigDist','normalize',1,'calchankelevery',1,'npts4h',0}; ...
%          {'sst3defaults',1,'dfunname','distAvgEigDist','normalize',1,'calchankelevery',1,'npts4h',-1}; ...
%          {'sst3defaults',1,'dfunname','distAvgEigDist','normalize',0,'calchankelevery',0,'npts4h',0}; ...
%          {'sst3defaults',1,'dfunname','distAvgEigDist','normalize',0,'calchankelevery',0,'npts4h',-1}; ...
%          {'sst3defaults',1,'dfunname','distAvgEigDist','normalize',0,'calchankelevery',1,'npts4h',0}; ...
%          {'sst3defaults',1,'dfunname','distAvgEigDist','normalize',0,'calchankelevery',1,'npts4h',-1}; ...

%          {'sst3defaults',1,'dfunname','distWeightedEigDist','normalize',1,'calchankelevery',0,'npts4h',0}; ...
%          {'sst3defaults',1,'dfunname','distWeightedEigDist','normalize',1,'calchankelevery',0,'npts4h',-1}; ...
%          {'sst3defaults',1,'dfunname','distWeightedEigDist','normalize',1,'calchankelevery',0,'npts4h',0,'symnorm',0}; ...
%          {'sst3defaults',1,'dfunname','distWeightedEigDist','normalize',1,'calchankelevery',1,'npts4h',0}; ...
%          {'sst3defaults',1,'dfunname','distWeightedEigDist','normalize',1,'calchankelevery',1,'npts4h',-1}; ...
%          
%          {'sst3defaults',1,'dfunname','distWeightedEigDist','normalize',0,'calchankelevery',0,'npts4h',-1}; ...
%          {'sst3defaults',1,'dfunname','distWeightedEigDist','normalize',0,'calchankelevery',1,'npts4h',0}; ...
%          {'sst3defaults',1,'dfunname','distWeightedEigDist','normalize',0,'calchankelevery',1,'npts4h',-1}; ...
        };
    chNamesBasic={'True';'Inverse';'Random';'Ones';'Zeros';'SST'};
    
    save(sprintf('./Data/params.mat'),'useControlAlgorithms','nSeriesBase', 'useBalancedStats', 'T', ...
    'noiseLevel','gNoiseSigmaInitial','pOutlier', 'minChangeDistance', 'pchange', ...
    'w', 'nwdiff', 'n', 'points', 'thresholds', 'findLocsFun', ...
    'firstW', 'statNames','rsstParams','chNamesBasic','uniformNoises','gNoiseSigmas','range','scale2range');
end

nNoisesNormal=length(gNoiseSigmas);
nNoisesUniform=length(uniformNoises);
nSeries=nSeriesBase*nNoisesNormal*nNoisesUniform;

nTh=length(thresholds);
nStats=length(statNames);
nP=numel(points);
x=cell(nSeries,1);
changes=cell(nSeries,1);
chagneLocs=cell(nSeries,1);
locsT=cell(nSeries,1);
t=cell(nSeries,1);
xTrue=cell(nSeries,1);

nrsstParams=size(rsstParams,1);
if useControlAlgorithms
    nAlgs=length(rsstParams)+length(chNamesBasic);
else
    nAlgs=length(rsstParams)+1; % one for sst
end
times=zeros(nSeries,nAlgs);
statVals=zeros(nSeries,nTh+1,nP,nStats,nAlgs);    
thBest=zeros(nSeries,nP,nStats,nAlgs);    
cpsps=zeros(nSeries,nAlgs,nP);
esrs=zeros(nSeries,nAlgs,nP);
desrs=zeros(nSeries,nAlgs,nP);
aesrs=zeros(nSeries,nAlgs,nP);
esrsl=zeros(nSeries,nAlgs,nP,nTh+1);
desrsl=zeros(nSeries,nAlgs,nP);
aesrsl=zeros(nSeries,nAlgs,nP);
esrse=zeros(nSeries,nAlgs,nP+1);
desrse=zeros(nSeries,nAlgs,nP+1);
aesrse=zeros(nSeries,nAlgs,nP+1);



% generate data
if LoadSavedSignal
   load(sprintf('./Data/signals.mat'),'x','locsT','t','xTrue');
else
    nxt=1;
    for series=1:nSeriesBase
        fprintf('Generating Series %d of %d\n',series,nSeriesBase);
        [xbase,locsTbase,tbase,xTruebase]=produceSingle(T,'pchange',pchange,'minChangeDistance',minChangeDistance,'noiseLevel',noiseLevel,'poutlier',pOutlier,'noiseSigma',gNoiseSigmaInitial,'range',range,'scale2range',scale2range);
        xbase=xbase';
        for u=1:nNoisesUniform
            for g=1:nNoisesNormal
                x{nxt}=xbase+uniformNoises(u).*rand(size(xbase))+gNoiseSigmas(g).*randn(size(xbase));
                locsT{nxt}=locsTbase; t{nxt}=tbase; xTrue{nxt}=xTruebase;
                nxt=nxt+1;
            end
        end        
        %dlmwrite('./Data/x{series}.dat',x{series},'\n');
        %dlmwrite('./Data/xTrue{series}.dat',xTrue{series},'\n');
        %dlmwrite('./Data/t{series}.dat',t{series},'\n');
        %dlmwrite('./Data/locsT{series}.dat',locsT{series},'\n');
        %save('./Data/signal.mat');    
    end
    save(sprintf('./Data/signals.mat'),'x','locsT','t','xTrue');
end    
% find the changes according to fixed control algorithms 
if LoadSavedChanges
    load(sprintf('./Data/changes.mat'),'changes','chNames','times');
else
    for series=1:nSeries                            
        if useControlAlgorithms
            chNames=chNamesBasic;
        else
            chNames={'SST'};
        end        
        ch=[];
        current=1;    

        if useControlAlgorithms
            fprintf('(%02d of %02d)\tCalculating base change points .... ',series,nSeries);
            tmp=tic;
            y=t{series};
            tmp=toc(tmp); times(series,current)=tmp;
            ch=[ch;y];

            current=current+1;

            tmp=tic;
            tLocs=findLocsTh((1-t{series}),0.1,0.5);
            y=zeros(size(t{series}));
            y(tLocs)=1;
            tmp=toc(tmp); times(series,current)=tmp;
            ch=[ch;y];

            current=current+1;

            tmp=tic;
            y=rand(size(t{series}));
            tmp=toc(tmp); times(series,current)=tmp;
            ch=[ch;y];

            current=current+1;

            tmp=tic;
            y=ones(size(t{series}));
            tmp=toc(tmp); times(series,current)=tmp;
            ch=[ch;y];

            current=current+1;

            tmp=tic;
            y=zeros(size(t{series}));
            tmp=toc(tmp); times(series,current)=tmp;
            ch=[ch;y];

            current=current+1;
            fprintf('Done\n');        
        end
    % find the changes according to different sst approaches
        fprintf('(%02d of %02d)\tCalculating SST change points .... ',series,nSeries);
        tmp=tic;
        y=sst(x{series},w,n,3,0,n);
        tmp=toc(tmp); times(series,current)=tmp;
        ch=[ch;y];    
        current=current+1;

%         tmp=tic;
%         y=sst3(x{series},w,n,-1,-1,-1,-1,-1);
%         tmp=toc(tmp); times(series,current)=tmp;
%         ch=[ch;y];    
%         current=current+1;
        
        fprintf('Done\n')
        fprintf('(%02d of %02d)\tCalculating RSST change points: \n',series,nSeries);
        for i=1:nrsstParams
            fprintf('(%02d of %02d)\t\tCalculating RSST change points (%d of %d) ... ',series,nSeries,i,nrsstParams);
            tmp=tic;
            [y]=cpd(x{series},w,n,-1,deal(squeeze(rsstParams{i})));
            tmp=toc(tmp); times(series,current)=tmp;
            ch=[ch;y];         
            switch(i)
                case 1
                    chNames{current,1}=sprintf('MZ',i); 
                case 2
                    chNames{current,1}=sprintf('RSST',i); 
                otherwise
                    chNames{current,1}=sprintf('SSA%02d',i); 
            end
            current=current+1;            
            fprintf('Done\n');
        end


    % normalize data
        times=(1000.*times)./(numel(x{series}));       % find time per point in ms
        ch(isnan(ch))=0;
        ch=abs(ch);            
        for current=1:nAlgs
            mxchg=max(ch(current,:));
            if(mxchg>eps)
                ch(current,:)=ch(current,:)./mxchg;
            end
        end
        changes{series}=ch;
    end
    save(sprintf('./Data/changes.mat'),'changes','chNames','times');
end        
if LoadSavedChangeLocs
    load(sprintf('./Data/changeLocs.mat'),'changeLocs');
else
    for series=1:nSeries
        ch=changes{series};
    % now calcualte localization results    
        fprintf('(%02d of %02d)\tLocalizing Changes: \n',series,nSeries);
        chLocs=cell(nAlgs,nTh);
        for i=1:nAlgs
            fprintf('(%02d of %02d)\t\tLocalizing Changes in algorithm %d of %d: ',series,nSeries,i,nAlgs);
            chLocs(i,:)=findLocsFun(ch(i,:),thresholds,firstW);
            fprintf('Done\n'); 
        end 
        changeLocs{series}=chLocs;
    end
    save(sprintf('./Data/changeLocs.mat'),'changeLocs');
end

if LoadSavedTraditionalMetrics
    load(sprintf('./Data/cpquality.mat'),'statVals','thBest');
else
    for series=1:nSeries                    
    % evaluate traditional metrics    
        fprintf('(%02d of %02d)\tEvaluating Traditional Metrics :\n',series,nSeries);    
        for th=1:nTh
            fprintf('(%02d of %02d)\t\tMetrics for threshold %d of %d: ',series,nSeries,th,nTh);
            for p=1:nP;
                [stats,timeStats]=cpquality({chLocs{:,th}},locsT{series},T,points(p),useBalancedStats);
                s=cell2mat(struct2cell(stats));
                statVals(series,th,p,:,:)=s(1:11,:);                        
            end
            fprintf('Done\n');
        end        

    % find maximum stats for every algorithm in every condition 
        fprintf('(%02d of %02d)\tFinding Best Threshold for each metric ... ',series,nSeries);
        for a=1:nAlgs                
            for p=1:nP
                for s=1:nStats
                    [statVals(series,nTh+1,p,s,a),nBest]=max(squeeze(statVals(series,1:nTh,p,s,a)));
                    thBest(series,p,s,a)=nBest;
                end
            end
        end
        fprintf('Done\n');
    end
    save(sprintf('./Data/cpquality.mat'),'statVals','thBest');
end
if LoadSavedCPSP
    load(sprintf('./Data/cpsp.mat'),'cpsps');
else
    for series=1:nSeries                    
        ch=changes{series};
    % evaluate new metrics
        fprintf('(%02d of %02d)\tEvaluating CPSP: \n',series,nSeries);
        for p=1:nP    
            fprintf('(%02d of %02d)\t\tCPSP for nPoints %d of %d: ',series,nSeries,p,nP);
            [cpsps(series,:,p)]=cpsp(ch',t{series},points(p));
            fprintf('Done\n');
        end    
    end
    save(sprintf('./Data/cpsp.mat'),'cpsps');
end
if LoadSavedESR
    load(sprintf('./Data/esrs.mat'),'esrs','desrs','aesrs');
else
    for series=1:nSeries                    
        fprintf('(%02d of %02d)\tEvaluating ESR: \n',series,nSeries);
        ch=changes{series};
        for current=1:nAlgs        
            fprintf('(%02d of %02d)\t\tESR for algorithm %d of %d: ',series,nSeries,current,nAlgs);
            for p=points
                esrs(series,current,p+1)=esr(ch(current,:),t{series},p);
            end
            desrs(series,current,:)=squeeze(esrs(series,current,:))-[0;squeeze(esrs(series,current,1:end-1))] ;
            for j=1:nP
                aesrs(series,current,j)=sum(squeeze(esrs(series,current,1:j)))./j;
            end
            fprintf('Done\n');
        end            
    end
    save(sprintf('./Data/esrs.mat'),'esrs','desrs','aesrs');
end
if calcESRLoc
if LoadSavedESRLoc
    load(sprintf('./Data/esrsloc.mat'),'esrsl','desrsl','aesrsl');
else
    for series=1:nSeries                    
        fprintf('(%02d of %02d)\tEvaluating ESR (after Localization): \n',series,nSeries);              
        for current=1:nAlgs        
            fprintf('(%02d of %02d)\t\tESR for algorithm %d of %d (After Localization): ',series,nSeries,current,nAlgs);
            for th=1:nTh
                ch=zeros(size(changes{series}(current,:))); ch(changeLocs{series}{current,th})=1;                
                for p=points
                    esrsl(series,current,p+1,th)=esr(ch,t{series},p);
                end                                
            end
            esrsl(series,current,:,nTh+1)=max(squeeze(esrsl(series,current,:,1:nTh)),[],2);
            fprintf('Done\n');
        end        
        for current=1:nAlgs
            fprintf('(%02d of %02d)\t\tESRD/L for algorithm %d of %d (After Localization) ... ',series,nSeries,current,nAlgs);
            desrsl(series,current,:)=squeeze(esrsl(series,current,:,nTh+1))-[0;squeeze(esrsl(series,current,1:end-1,nTh+1))] ;
            for j=1:nP
                aesrsl(series,current,j)=sum(squeeze(esrsl(series,current,1:j,nTh+1)))./j;
            end
            fprintf('Done\n');
        end
    end
    save(sprintf('./Data/esrsloc.mat'),'esrsl','desrsl','aesrsl');
end
end
if calcESRN
if LoadSavedESRN
    load(sprintf('./Data/esrse.mat'),'esrse','desrse','aesrse');
else
    for series=1:nSeries                    
        fprintf('(%02d of %02d)\tEvaluating ESR (Exponential Scaling): \n',series,nSeries);
        ch=changes{series};
        for current=1:nAlgs        
            fprintf('(%02d of %02d)\t\tESR (Exponential) for algorithm %d of %d: \n',series,nSeries,current,nAlgs);
            for p=1:nP
                fprintf('(%02d of %02d)\t\t\tESR (Exponential) for algorithm %d of %d (nPoints %d of %d): ',series,nSeries,current,nAlgs,p,nP+1);
                esrse(series,current,p)=esr(ch(current,:),t{series},points(p),'normalScaling',points(end));
                fprintf('Done\n');
            end
            fprintf('(%02d of %02d)\t\t\tESR (Exponential) for algorithm %d of %d (nPoints %d of %d<INF>): ',series,nSeries,current,nAlgs,nP+1,nP+1);
            esrse(series,current,nP+1)=esr(ch(current,:),t{series},inf,'expScaling');
            fprintf('Done\n');
            desrse(series,current,:)=squeeze(esrse(series,current,:))-[0;squeeze(esrse(series,current,1:end-1))] ;
            for j=1:nP+1
                aesrse(series,current,j)=sum(squeeze(esrse(series,current,1:j)))./j;
            end            
        end            
    end
    save(sprintf('./Data/esrse.mat'),'esrse','desrse','aesrse');
end       
end

if displayResults    
% display things
    dispTestResultsCPD
end
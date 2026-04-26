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

%rng(pi);
clear; 
warning('off','stats:regress:RankDefDesignMat');
load('./Data/ccausalitysynthetic.mat');
%nModels=min(size(Ms,1),2);
nSteps=100;         %threshold steps

maxLag=delayRange(2)+1;
w=min(5,maxLag);
w_n=w+1;

detector=   {'None','Optimal','RSST'};
method=     {'Granger','DCD'};
cfroml=[0,0,1];
thresholdRange=[0,1; ...
                0,1; ...
                0,1 ...
                ];
%needsMaxMinDet=[1;0;1];
%accept0And1=[1;0;1];
applicableDet=[1,1,1; ...
               0,1,1; ...               
              ];
fields={'F1','MCC'};

nm=length(method);
nd=length(detector);
nF=length(fields);



Mg=cell(nm,nd,nModels);
c=cell(nModels,nd);
crsst=cell(nModels,1);
lrsst=cell(nModels,1);
locs=cell(nModels,nd);
sg=cell(nm,nd,nModels);
st=cell(nm,nd,nModels);
rocg=cell(nm,nd,nModels);
statF=cell(nm,nd,nModels,nF);
MgSave=cell(nm,nd,nModels,nF);
thSaved=cell(nm,nd,nModels,nF);
bestStat=cell(nm,nd,nModels,nF);
mytimes=zeros(nm,nd,nModels);
detTimes=zeros(nModels,nd);
rdet=0;
for model=1:nModels
    fprintf('Model number %d\n',model);
    M=Ms{model}; E=Es{model}; X=xs{model}; ndim=ndims{model}; chgLoc=chgLocs{model}; chgs=chgss{model};    
    nx=size(X,1);
    fprintf('\tCalculating change points ');
    for det=1:nd
        fprintf('%s',(detector{det}));
        switch(lower(detector{det}))
            case 'none'
                startT=tic;
                c{model,det}=X;
                detTimes(model,det)=toc(startT);
            case 'optimal'
                startT=tic;
                c{model,det}=chgs;
                locs{model,det}=chgLoc;
                detTimes(model,det)=toc(startT);
            case 'rsst'
                c{model,det}=zeros(nx,size(X,2));
                if(cfroml(det))
                    if isempty(lrsst{model})
                        startT=tic;
                        for i=1:nx
                            [crsst{model}(i,:),locs{model,det}{i}]=rsst(X(i,:),w,w_n,-1,'th',noiseSigma+1e-12);                            
                            c{model,det}(i,:)=addGaussiansInLocs(locs{model,det}(i),size(X,2),maxLag/3);
                        end
                        detTimes(model,det)=toc(startT);
                        lrsst{model}=locs{model,det};
                        rdet=det;
                    else
                        startT=tic;
                        for i=1:nx                            
                            c{model,det}(i,:)=addGaussiansInLocs(lrsst{model}(i),size(X,2),maxLag/3);
                        end
                        detTimes(model,det)=toc(startT)+detTimes(model,rdet);
                    end
                else
                    if isempty(crsst{model})
                        startT=tic;
                        for i=1:nx
                            [c{model,det}(i,:),locs{model,det}{i}]=rsst(X(i,:),w,w_n,-1,'th',noiseSigma+1e-12);
                        end
                        detTimes(model,det)=toc(startT);
                        lrsst{model}=locs{model,det};
                        rdet=det;
                    else
                        startT=tic;
                        c{model,det}=crsst{model};
                        detTimes(model,det)=toc(startT)+detTimes(model,rdet);
                    end
                end
        end
        fprintf('(%4.3fs) ... ',detTimes(model,det));
    end
    fprintf('...Done\n');
    for meth=1:nm
        for det=1:nd
            if(~applicableDet(meth,det))
                continue;
            end            
            fprintf('\t%s (%s) ...',method{meth},detector{det});            
%             if((thresholdRange(meth,2)-thresholdRange(meth,1))<1/nSteps)
%                 switch(lower(detector{det}))
%                     case 'none'
%                         [dummy,F]=detectCC(X,0,'method',method{meth},'cdmethod',detector{det},'maxlag',maxLag,'c',X);
%                     case 'optimal'
%                         [dummy,F]=detectCC(X,0,'method',method{meth},'cdmethod',detector{det},'maxlag',maxLag,'c',chgs,'locs',chgLoc);
%                     case 'rsst2'
%                         [dummy,F]=detectCC(X,th,'method',method{meth},'cdmethod','rsst','maxlag',maxLag,'cfroml',1);
%                     otherwise                                            
%                         [dummy,F]=detectCC(X,0,'method',method{meth},'cdmethod',detector{det},'maxlag',maxLag);
%                 end
%                 mxF=max(max(F));   mnF=min(min(F));                                                
%             else
%                 mxF=thresholdRange(meth,2);  mnF=thresholdRange(meth,1);
%             end
            mxF=thresholdRange(meth,2);  mnF=thresholdRange(meth,1);
            step=(mxF-mnF)/nSteps;
            th=mnF:step:mxF;
            startTic=tic;
            [Mg{meth,det,model},F,cvUsed]=detectCC(X,th,'method',method{meth},'cdmethod','none','maxlag',maxLag,'c',c{model,det},'locs',locs{model,det});
%             switch(lower(detector{det}))
%                 case 'none'
%                     [Mg{meth,det,model},F,cvUsed]=detectCC(X,th,'method',method{meth},'cdmethod',detector{det},'maxlag',maxLag,'c',X);
%                 case 'optimal'
%                     [Mg{meth,det,model},F,cvUsed]=detectCC(X,th,'method',method{meth},'cdmethod',detector{det},'maxlag',maxLag,'c',chgs,'locs',chgLoc);
%                 case 'rsst2'
%                     [Mg{meth,det,model},F,cvUsed]=detectCC(X,th,'method',method{meth},'cdmethod','rsst','maxlag',maxLag,'cfroml',1);
%                 otherwise                                            
%                     [Mg{meth,det,model},F,cvUsed]=detectCC(X,th,'method',method{meth},'cdmethod',detector{det},'maxlag',maxLag);
%             end
            totalTime=(toc(startTic)+detTimes(model,det))/(T*numel(th));
            mytimes(meth,det,model)=totalTime;
            rocg{meth,det,model}=zeros(numel(th),3);
            statF{meth,det,model}=zeros(numel(th),1);            
            for i=1:numel(th)                                    
                [sg{meth,det,model},st{meth,det,model}]=cmpModels(M,Mg{meth,det,model}(:,i));
                rocg{meth,det,model}(i,:)=[1-sg{meth,det,model}.specificity,sg{meth,det,model}.sensitivity,th(i)];                                        
                for fld=1:nF                    
                    fldVal=getfield(sg{meth,det,model},lower(fields{fld}));
                    statF{meth,det,model,fld}(i)=fldVal;                    
                    if(th(i)<eps); continue; end;
                    if(th(i)>1-eps); continue; end;
                    if isempty(bestStat{meth,det,model,fld}) || fldVal>bestStat{meth,det,model,fld}                        
                        MgSave{meth,det,model,fld}=Mg{meth,det,model}(:,i);
                        thSaved{meth,det,model,fld}=th(i);                
                        bestStat{meth,det,model,fld}=fldVal;                        
                    end                    
                end          
            end
            fprintf('... DONE in %4.3f ms/point(%8.3fHz)\n',1000*mytimes(meth,det,model),1/mytimes(meth,det,model));
%            fileName=sprintf('%s_%s_%03d',meth,det,model);
       end
    end
end

dispAndSaveSynthetic
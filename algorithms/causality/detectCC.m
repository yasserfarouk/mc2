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

function [M,s,cv,c]=detectCC(x,th,varargin)
% reads a multidimensional time series and detects the causality structure
% of changes in them .
%
% Inputs:
% ======
% x             n*T array of time series to find the relations between its
%               n timeseries
% th            the threshold used to get the granger causality. can be set
%               to a vector of length m and in this case the output M will be
%               a cell array of n*m members where M(:,i) corresponds to
%               th(i)
% varargin      see the swich statement
%
%
% Outputs:
% =======
% M             The model learned
% s             The value of the statisitc used in the detection (if any)
% cv            The value of the critical value overwhich s was considered significant and a relation is created    
% c             The change scores of the input signal (x) 
% 
%
    findSelfLoops=boolean(0);
    s=[]; cv=[];
    method='granger';
    cdmethod='rsst';
    cdsparams=[];
    cdoparams=[];
    c=[];
    locs=[];
    cFroml=0;
    lFromc=0;
    maxLag=10;
    locTh=0.5;
    locFun=@findLocsThMean;
    locFunParams=[];
    locSigma=maxLag/2;
    nArgs=size(varargin,2);
    if(nArgs>0)
        if(mod(nArgs,2)~=0)
            error('The optional arguments must be in the form name,value so they must be even!!!');
        end
        for i=1:2:nArgs
            switch(lower(varargin{i}))                   
                case {'locsigma','localizationstddev'}        % the standard deviation of the gaussians to be added around every change point to construct the change score (default maxLag/3)
                    locSigma=varargin{i+1};
                case {'cdalg','chagnedetectoralgorithm','cdmethod'} % the change detection algorithm to be used. default is 'rsst'
                    cdmethod=varargin{i+1};
                case {'cdparams','cdstdparams'}             % a cellarray of n*1 parameters to be passed to the change detector in the form of standard arguments
                    cdsparams=varargin{i+1};
                case {'cdoptionalparams'}                   % a cell array of n*2 (name,value) to be passed to the change detector in the form of varargin
                    cdoparams=varargin{i+1};
                case {'maxlag','maxdelay'}              % The maximum delay to look for relations in the time series dynamics (default 10)
                    maxLag=(varargin{i+1});                        
                case {'findselfloops','fsl'} 
                    findSelfLoops=(varargin{i+1});                        
                case {'algorithm','method','alg'}   % the algorithm used for detecting causality out of the signal
                    method=(varargin{i+1});       
                case {'changescores','cs','c'}      % a matrix of same dimensions as x giving the change points to be used. default found form rsst
                    c=(varargin{i+1});              
                case {'changelocations','locs'}     % a cell array with one item for each time series containing the locations of change points. default is to be found from rsst
                    locs=(varargin{i+1});              
                case {'cfroml','changefromlocs'}    % if nonzero then c will be ignored and created based on Locs. default 0
                    cFroml=(varargin{i+1});
                case {'lfromc','locsfromchange'}    % if nonzero then Locs will be ignored and created based on c. default 0
                    lFromc=(varargin{i+1});
                case {'locth'}                      % the threshold to be used to localize the changes ( useful only if lfromc was nonzero and locs was empty) default is 0.5
                    locTh=(varargin{i+1});
                case {'locfun'}                     % useful only if lfromc was nonzero and locs was empty. this function receives the result and M and K and may receive other parameters if pFunParams was given
                    locFun=varargin{i+1};
                case {'locfparams','locparams','locfunparams'} % any parameters to be passed to locfun
                    locFunParams=varargin{i+1};
                otherwise
                    error('Unknown argument: %s',varargin{i});
            end
        end
    end    
    method=lower(method);
    cdmethod=lower(cdmethod);
    locSigma=maxLag/2;
    % initialize the change detector's standard parameters if not given
    if(strcmp(cdmethod,'rsst') && isempty(cdsparams))
        cdsparams={maxLag;maxLag+1;-1};
    end
    
    % detect c and locs as needed
    if(cFroml)
        [c,locs]=changeDetect(x,[],locs,0,needsLoc(method) || needsC(method),cdmethod,cdsparams,cdoparams);
    else
        if(lFromc)
            [c,locs]=changeDetect(x,c,[],needsLoc(method) || needsC(method),0,cdmethod,cdsparams,cdoparams);
        else
            [c,locs]=changeDetect(x,c,locs,needsC(method),needsLoc(method),cdmethod,cdsparams,cdoparams);
        end
    end

    n=size(x,1);
    % confirm that you have locs if you need it
    if(needsLoc(method) && isempty(locs))
        if isempty(c)
            error('Locs is to be found from c but c is not available');
        end    
        for i=1:n
            locs{i}=locFun(c(i,:),locTh,locFunParams);
        end
    end

    % confirm that you have c if you need it
    if(needsC(method) && isempty(c))
        if isempty(locs)
            error('c is to be found from locs but locs is not available');
        end
        c=addGaussiansInLocs(locs,size(x,2),locSigma);
    end


    switch(lower(method))
        case 'granger'
            %[M,s,cv]=detectGC(c,0,maxLag,'findSelfLoops',findSelfLoops,'th',th);            
            [M,s,cv]=detectGC(c,th,maxLag,'findSelfLoops',findSelfLoops);            
        case 'normality'
            [M,s]=detectCCUsingNormalityTest(locs,th,'fsl',findSelfLoops,'maxLag',maxLag);
        case 'dcd'      % stands for delay consistency discrete
            [M,s]=detectCCUsingDelayConsistency(locs,th,'fsl',findSelfLoops,'maxLag',maxLag);
        otherwise
            error('unknown method(%s)',method);
    end
end

function n=needsC(method)    
    switch(lower(method))
        case 'granger'
            n=1;
        case 'normality'
            n=0;
        case 'dcd'
            n=0;
        otherwise
            error('unknown method(%s)',method);
    end
end
function n=needsLoc(method)    
    switch(lower(method))
        case 'granger'
            n=0;
        case 'normality'
            n=1;
        case 'dcd'
            n=1;
        otherwise
            error('unknown method(%s)',method);
    end
end


function [c,locs]=changeDetect(x,c,locs,needC,needL,cdmethod,cdsparams,cdoparams)
    n=size(x,1);    
    
    % correct needs
    needC=needC && (isempty(c));
    needL=needL && (isempty(locs));
    
    if(needL);     locs=cell(n,1);      end;
    if(needC);     c=zeros(size(x));    end;
    if(~needL && ~needC)
        return;
    end
    % decide how many outputs to get from the detector
    switch (cdmethod)
        case 'rsst'
            needF=needC; needS=needL;
        otherwise
            error('Unknown change point detector: %s',cdmethod);
    end
    
    for i=1:n
        if needC
            % call detector
            [first,second]=callDetector(x(i,:),cdmethod,cdsparams,cdoparams,needF,needS);    
            %decide which output is which
            switch (cdmethod)
                case 'rsst'
                    if(needC);    c(i,:)=first;  end;
                    if(needL);    locs{i}=second; end;
                otherwise
                    error('Unknown change point detector: %s',cdmethod);
            end
        end
        if needL
            locs{i}=findLocsTh(c(i,:));
        end
    end
    
end

function [c,l]=callDetector(x,cdmethod,cdsparams,cdoparams,needF,needS)    
    c=[];l=[];dummy=[];
    if(needF && needS)
        callstr=sprintf('[c,l]=%s(x',cdmethod);    
    elseif needF
        callstr=sprintf('[c]=%s(x',cdmethod);
    elseif needS
        callstr=sprintf('[dummy,l]=%s(x',cdmethod);
    else
        callstr=sprintf('%s(x',cdmethod);
    end    
    if ~isempty(cdsparams) || ~isempty(cdoparams)
        callstr=[callstr,','];
    end
    if ~isempty(cdsparams)
        n=size(cdsparams,1);
        for i=1:n
            if isempty(cdsparams{i})
                callstr=[callstr,'[]'];
            elseif ischar(cdsparams{i})
                callstr=[callstr,'''',cdsparams{i},''''];
            else
                callstr=[callstr,'[', num2str(cdsparams{i}),']'];
            end
            if(i<n)
                callstr=strcat(callstr,',');
            end
        end
    end
    if ~isempty(cdoparams)
        n=size(cdoparams,1);
        for i=1:n
            if isempty(cdoparams{i,2})
                callstr=[callstr,'''',cdoparams{i,1},'''',',[]'];
            elseif ischar(cdsparams{i})
                callstr=[callstr,'''',cdoparams{i,1},'''',',''',cdoparams{i,2},''''];
            else
                callstr=[callstr,'''',cdoparams{i,1},'''',',[',num2str(cdoparams{i,2}),']'];
            end
            if(i<n)
                callstr=strcat(callstr,',');
            end
        end
    end    
    callstr=strcat(callstr,');');
    eval(callstr);
end

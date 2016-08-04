function locs=findLocsThin(r,varargin)
% first parameter is the time series in which we want to localize local
% maxima
% optional parameters in order
% th        The threshold underwhich we do not consider any thing (default
%           is 1e-12). If multiple thresholds are given then the result is
%           a cell array of changes with each threshold alone
% firstW    if firstW<-1 then the location of the local maxima is just returned
%           if 0>firstW>=-1 then one of the points to attain the local maxima is
%           to be used as the location of the maxima using:
%           firstW*firstLocalMaxima+(1-firstW)*lastLocalMaxima
%           where firstLocalMaxima(lastLocalMaxima) are the location within
%           a local vecinity of a local maxima that are within th value
%           from it
%           if positive then the location of the maxima is calculated as:
%           firstW*first+(1-firstW)*last 
%           where first and last are the points at which the signal goes
%           over and returns under the threshold th
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


    th=1e-6;    
    firstW=-2;     % default is location of maximum change
    nArgs=size(varargin,2);
    if(nArgs>0)
        th=varargin{1};
    end        
    if(nArgs>1)
        firstW=varargin{2};        
    end 
    m=length(th);
    if m==1
        thinned=ones(size(r)); thinned(r<th)=0;
        locs=findLocsTh(thinned,0.5,firstW);
    else
        locs=cell(1,m);
        for t=1:m
            thinned=ones(size(r)); thinned(r<th(m))=0;            
            locs{t}=findLocsTh(thinned,0.5,firstW);
        end        
    end    
end

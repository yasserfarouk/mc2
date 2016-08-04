function [r,u,s]=tspca(x,varargin)
% converts an M dimensional time series into a 1 dimensional time series by
% applying PCA to it. We assume that each dimension is along a row which
% means that each column corresponds to M values for a single time step and
% the matrix is M*T where T is the length of the time series.
%
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%

normalize=1;
doPlot=0;
M=size(x,1);
T=size(x,2);
rotated=false;
if M>T
    rotated=true;
    x=x';
    M=size(x,1);
    T=size(x,2);
end
maxT=[];
nArgs=size(varargin,2);
if(nArgs>0)
    if(mod(nArgs,2)~=0)
        error('The optional arguments must be in the form name,value so they must be even!!!');
    end
    for i=1:2:nArgs
        switch(lower(varargin{i}))       
            case {'maxt'}                   % the maximum number of time steps to use in the SVD. These will be taken at equal distances from 1:T. (default max(20*M,log T))
                maxT=(varargin{i+1});
            case {'plot'}                   % whether or not to plot the signal after normalization
                doPlot=(varargin{i+1});                           
            case {'normalize','norm'}       % if nonzero each dimension will be z normalized by subtracting the mean and dividing by the std. dev.
                normalize=(varargin{i+1});               
            otherwise
                error('unknown command: %s',varargin{i});
        end
    end
end
if isempty(maxT)
    maxT=max(100*M,ceil(log2(T)));
end
if(normalize)
    for i=1:M
        m=mean(x(i,:));
        s=std(x(i,:));
        if abs(s)>1e-6*m
            x(i,:)=(x(i,:)-m)./s;
        else
            x(i,:)=x(i,:)-m;
        end
    end
end
if doPlot
    plot(x');
end
step=max(1,ceil(T/maxT));
if step==1
    [u,s,r]=svd(x(:,1:step:T),'econ');
    s=diag(s);
    r=r(:,1)';
else
[u,s,v]=svd(x(:,1:step:T),'econ');
s=diag(s);
r=u(:,1)'*x;
if rotated
    r=r';
end
end

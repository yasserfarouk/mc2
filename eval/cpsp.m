function q=cpsp(x,t,nPoints)
% finds the quality of change points in x using t as true change point
% locaitons. Both x and t must be the same size
%
% x         output of change point algorithm. If x is n*m matrix then it is
%           treated as m outputs of m algorithms with a time series of
%           length n and nPoints must be n points. This means the each
%           algorithm's output must be in a column not a row
% t         true change points (e.g. 1 at change points and 0 otherwise)
%           the value at each change point may be treated as weight and the
%           more this value is the more important this change point will be
%           in the result. This must be a vector
% nPoints   width around change points that is assumed OK. If a negative
%           number it is selected as distance to next true chagne point
%           divided by -nPoints
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



% set constants
epsilon=1e-6;
% if(nargin<4)
%     epsilon=1e-6;
% else
%     epsilon=th;
% end

% check input
if(nargin<3)
    nPoints=0;
end

if(size(x,1)==1) && (size(x,2)~=1)
    x=x';
end
if ~isvector(t)
    error('t must be a vector');
end
m=size(x,2);
n=size(x,1);
if n ~= numel(t)
    error('The number of elements in t must be equal to the number of columns in x');
end

% normalize
t=t./max(t);
q=zeros(m,1);
locs=find(abs(t)>epsilon);
p=length(locs);

y=max(x);
for i=1:m
    x(:,i)=x(:,i)./y(i);
end

for c=1:m       % for each column of x        
    q(c)=mychangeQuality(x(:,c),n,locs,p,nPoints);    
end


end

% x here MUST be a vector and assumes that its sum is 1
function q=mychangeQuality(x,xn,locs,p,nPoints)
q=0;
% evaluate    
counted=zeros(size(x));
for i=1:p
    % find nearest change point before and after taking care of
    % boundareies
    if(i==1)
        nearestB=1;        
    else
        nearestB=locs(i-1);
    end
    if (i==p) 
        nearestA=xn;
    else
        nearestA=locs(i+1);
    end
    % set the number of points to count around the t at each direction
    if nPoints>=0
        nPB=nPoints;
        nPA=nPoints;
    else
        nPA=round((nearestA-locs(i))/nPoints);
        nPB=round((locs(i)-nearestB)/nPoints);
    end
    % find the range of points to count
    range=[max(1,locs(i)-nPB) min(xn,locs(i)+nPA)];

    % do count
    for j=range(1):range(2)
        if ~counted(j)
            q=q+x(j);
            counted(j)=1;
        end
    end                
end    
% normalize q
if(sum(x)>eps)
    q=q./sum(x);
else
    q=0;
end
end
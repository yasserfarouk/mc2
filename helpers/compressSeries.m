function s=compressSeries(d,n,pointsAreRows)%,mx,mn)
% changes the series length to be equal to n. can be applied to vectors
% or 2d matrices
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

if isvector(d)
    m=numel(d);
    if m ==n
        s=d;
    else
        if(size(d,1)==1)
            s=zeros(1,n);
        else
            s=zeros(n,1);
        end
        for i=1:n
            s(i)=d(max(1,floor(i*m/n)));
        end            
    end
else
    if nargin<3
        pointsAreRows=1;    
    end
    m=size(d,pointsAreRows+1);
    if m ==n
        s=d;
    else    
        if pointsAreRows
            s=zeros(size(d,1),n);
            for i=1:n
                s(:,i)=d(:,max(1,floor(i*m/n)));
            end    
        else
            s=zeros(n,size(d,2));
            for i=1:n
                s(i,:)=d(max(1,floor(i*m/n)),:);
            end    
        end
    end 
end
%s=((mx-mn)./(max(s)-min(s))).*s+mn;
end

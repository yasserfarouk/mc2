function [V,E,kind]=drawModel(M,n)
% draws a causality model using grPlot (M is the model and n (optional) the
% number of processes)
showDelay=0;
nSeries=size(M,1);
if(nargin<2); n=0; end;
n=max(nSeries,n);    
V=[cos(2*pi*[1:n]'/n) sin(2*pi*[1:n]'/n)]; 
kind='d';

E=[];
for i=1:nSeries    
    nI=size(M{i,1});
    for j=1:nI                
        if(showDelay)
            E=[E;M{i,1}(j,1),i,round(M{i,1}(j,2))];        
        else
            E=[E;M{i,1}(j,1),i];        
        end
    end
end
grPlot(V,E,kind,'%d','%d');
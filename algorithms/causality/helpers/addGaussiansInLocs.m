function x=addGaussiansInLocs(locs,T,locSigma,weights)
% generates a matrix n*T of real values where n is the length(locs) and
% each one in locs is an array of locatiosn to put the gaussians in.
% weights is exactly the same as locs but gives the weight of the gaussian
% at every point. if no weights are given then all gaussians will have the
% same weight
n=length(locs);
x=zeros(n,T);
if(nargin<4)
    weights=[];
end
nor=normpdf(-T:T,0,locSigma);
for i=1:n
    if(~isempty(locs{i}))
        m=numel(locs{i});
        T=size(x,2);        
        if isempty(weights)
            for j=1:m
                x(i,:)=x(i,:)+nor(T-locs{i}(j)+2:2*T-locs{i}(j)+1);
            end
        else
            for j=1:m
                x(i,:)=x(i,:)+weights{i}(j)*normpdf(1:T,locs{i}(j),locSigma);
            end
        end
    end
end
function [x,locs,y]= embedInside(x,locs,y,L,randRangeBefore,gNoiseSigma,onlyAfter,useSelfForExtension)
%EMBEDINSIDE Embedds a timeseries into another with specific length
%   See generatePatterns() for the meaning of all parameters
%
% onlyAfter     if true then x will be embedded in the beginning so that locs
%               does not change at all
% useSelfForExtension   if nonzero rather than extending using random data.
%               We use data from the same timeseries x

if~exist('randRangeBefore','var')
    randRangeBefore=[-0.5 0.5];
end
if~exist('gNoiseSigma','var')
    gNoiseSigma=0;
end
if~exist('onlyAfter','var')
    onlyAfter=false;
end

if~exist('useSelfForExtension','var')
    useSelfForExtension=true;
end


nAll=numel(x);
if L<nAll
    error('Input timeseries is longer then required output!!!');
end
nMotifs=numel(locs);
occurences=zeros(nMotifs,1);
for i=1:nMotifs
    occurences(i)=size(locs{i},1);
end
fl=locs{1}(1,2)-locs{1}(1,1);
smoothingFactor=min(floor(fl/8),max(5,floor(fl/15)));

if nAll<L
    if onlyAfter
        before=0;
    else
        before=floor((L-nAll).*rand(1,1));
    end
   after=L-nAll-before;
   if before>0
        if useSelfForExtension
            ind=randi(nAll,1,before);
            z=x(ind);
        else
            z=(randRangeBefore(2)-randRangeBefore(1))*rand(1,before)+randRangeBefore(1);             
        end
        z=smooth(z,smoothingFactor)';
        y=[zeros(size(z)) y];
        z=z+gNoiseSigma*randn(size(z));
        x=[z x];
        for i=1:nMotifs
            nz=numel(z);
            locs{i}=locs{i}+nz*ones(occurences(i),2);                
        end
   end
   if after>0
       if useSelfForExtension
           ind=randi(nAll,1,after);
           z=x(ind);
       else
           z=(randRangeBefore(2)-randRangeBefore(1))*rand(1,after)+randRangeBefore(1); 
           z=smooth(z,smoothingFactor)';
       end
       y=[y zeros(size(z))];
       z=z+gNoiseSigma*randn(size(z));
       x=[x z];
   end
end

end


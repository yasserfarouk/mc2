function [ Mc ] = correctCausalityStructure( M )
%corrects the causality structure by removing links that do not represent
%genuine direct causality relations
%
% Inputs:
% ======
% M         Causality structure to be optimized. See produceMulti for its
%           internal structure. Only the first dimension of the cell is
%           used. the rest are used for generation only
%
%
%
minDelay=4;
maxDelay=20;
Mc=M;
n=size(M,1);
for i=1:n
    m=size(Mc{i,1},1);
    toremove=[];
    for j=1:m
        p=Mc{i,1}(j,1);
        if M{i,1}(j,2)<minDelay || M{i,1}(j,2)>maxDelay
            toremove=[toremove; j];
        end        
    end
    Mc{i,1}(toremove,:)=[];
end

end


function [ locs,dists,executionTime,means,stats ] = parseMKoutput( str,x,lengthFirst )
%Parses the output of any exact motif discovery algorithm in the toolbox
%       
    if ~exist('lengthFirst','var')
        lengthFirst=false;
    end
    A=str2num(str);
    executionTime=A(end,end);
    A(end,:)=[];
    n=size(A,1);
    locs=cell(n,1);    
    if(nargout>2)
        means=cell(n,1);
        stats=cell(n,1);
    end
    dists=A(:,4);    
    if lengthFirst
        A=A(:,[2,3,1]);
    else
        A=A(:,1:3);
    end
    
    A=int32(A);
    for i=1:n
        locs{i}=[A(i,1)+1,A(i,1)+A(i,3);...
                 A(i,2)+1,A(i,2)+A(i,3)];
       if(nargout>3)
            [stats{i},means{i}]=calcStatsAndMean(locs{i},x,@abs);
       end
    end   
    
end


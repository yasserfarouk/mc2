function [ locs,maxSmallDistance,means,stats ] = parseMOENoutput( str,x )
%Parses the output of any exact motif discovery algorithm in the toolbox
%   
    A=str2num(str);    
    n=size(A,1);
    locs=cell(n,1);
    if(nargout>2)
        means=cell(n,1);
        stats=cell(n,1);
    end
    maxSmallDistance=max(A(:,4));
    A=A(:,1:3);
    A=int32(A);
    for i=1:n
        locs{i}=[A(i,2),A(i,1)+A(i,2)-1;...
                 A(i,3),A(i,1)+A(i,3)-1];
       if(nargout>1)
            [stats{i},means{i}]=calcStatsAndMean(locs{i},x,@abs);
       end
    end   
    
end


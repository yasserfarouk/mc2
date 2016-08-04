% This file is a part of the MC2 toolbox developed by Y. Mohammand and T. Nishida.
%Please do not remove this comment
%
% Using this file is governed by the license of MC2 which you can find in LICENSE.md
% 
% You can find more information about this toolbox here:
% - Yasser Mohammad and Toyoaki Nishida, "MC2: An Integrated Toolbox for Change, Causality, 
%   and Motif Discovery", 29th International Conference on Industrial, Engineering & 
%   Other Applications of Applied Intelligent Systems (IEA/AIE) 2016, pp. 128 -- 141.
% - Yasser Mohammad and Toyoaki Nishida, "Data Mining for Social Robotics", Springer 2016.
%

function [stop,toreuse]=stopIfAboveMean(x,begends,newbegends,toreuse,optional)
% Stops motif extension (internal function)
%
%
% For more information please consult the following publications: 
% ===============================================================
% Yasser Mohammad and Toyoaki Nishida, G-SteX: Greedy Stem Extension for
% Free-Length Constrained Motif Discovery, IEA/AIE 2012
% Yasser Mohammad and Toyoaki Nishida, CPMD: A Matlab Toolbox for Change
% Point and Constrained Motif Discovery, IEA/AIE 2012 
%
% Please cite the above mentioned publications if you are using this
% routine for your research.
%
    T=length(x);
    if max(max(newbegends))>T
        stop=1;
        return;
    end
    significance=0.05;
    epsilon=1e-20;    
    dPointDiffFun=@abs;
    acceptableMeanIncrease=0.1;
    if nargin>4 && ~isempty(optional)
        significance=optional{1};
        lop=length(optional);
        if lop>1
            epsilon=optional{2};
        end
        if lop>2
            dPointDiffFun=optional{3};
        end
        if(lop>3)
            acceptableMeanIncrease=optional{4};
        end
    end
    if ~isempty(begends)
        n=size(begends,1);
        m=max(begends(:,2)-begends(:,1)+1);
        X=zeros(n,m);       % convert to a cell if the sizes may not be equal
        for j=1:n
            X(j,:)=(x(begends(j,1):begends(j,2)));
        end
        % @todo here may be we can try to slide the ocurrences to get minimum
        % mean distance overall before calculating the means and variances.
        % his will requuire reading dFun as an input. The same can be done
        % in the other loop calculating newX!!
        d=0; d2=0; l=0;
        for i=1:n-1
            for j=i+1:n                
                diff=dPointDiffFun(X(i,:)-X(j,:));
                l=l+length(diff);
                d=d+sum(diff);
                d2=d2+sum(diff.*diff);
            end
        end
        meandist=d/l;
        vardist=(d2-(d^2/l))/(l-1);        
    else
        if isempty(toreuse)
        else
            meandist=toreuse{1}; vardist=toreuse{2}; l=toreuse{3}; X=toreuse{4};        
            begends=toreuse{5};        
            n=size(begends,1);        
        end
    end
    
    oldn=n;
    n=size(newbegends,1);
    m=max(newbegends(:,2)-newbegends(:,1)+1);
    newX=zeros(n,m);       % convert to a cell if the sizes may not be equal
    for j=1:n
        newX(j,:)=(x(newbegends(j,1):newbegends(j,2)));        
    end
    diffs=[]; newl=0;
    d=0; d2=0;
    new=cell(n,1);
    for i=1:oldn
        new{i}=setdiff(newbegends(i,1):newbegends(i,2),...
            begends(i,1):begends(i,2))-(newbegends(i,1)-1); 
        % [newbegends(i,1):begends(i,1)-1, begends(i,2)+1:newbegends(i,2)]-(newbegends(i,1)-1);
        new{i}(new{i}<1)=[];
        new{i}(new{i}>T)=[];
    end
    for i=1:n-1
        for j=i+1:n                
            diff=dPointDiffFun(newX(i,:)-newX(j,:));
            % add the new differences only to a list
            if i>oldn || j>oldn
                diffs=[diffs,diff];            
            else
                if(~isempty(new{i}) || ~isempty(new{j}))
                    diffs=[diffs,diff(union(new{i},new{j}))];
                end
            end
            newl=newl+length(diff);
            d=d+sum(diff);
            d2=d2+sum(diff.*diff);
        end
    end
    newmeandist=d/newl;
    newvardist=(d2-(d^2/newl))/(newl-1);
    %newvardist=((newmeandist^2)*newl+sum(d2-2.*newmeandist*d))/(newl-1);
    
    if newmeandist<=meandist*(1+acceptableMeanIncrease)
        stop=0;
    else
        if abs(vardist)<epsilon
            stop=(newmeandist-meandist>m*epsilon);
        else
            h=ztest(diffs,meandist,vardist,significance,'right');
            if isnan(h)
                warning('DGR:ZTESTNAN','z-test failed. we will just compare means');
                stop=(newmeandist-meandist>m*epsilon);
            else
                stop=h;
            end            
        end
    end
    
    if(stop)
        toreuse={meandist,vardist,l,X,begends};
    else
        toreuse={newmeandist,newvardist,newl,newX,newbegends};
    end
end

function locs=findLocsThMax(r,varargin)
    K=numel(r);
    locs=[];
    th=1e-12;
    nArgs=size(varargin,2);
    if(nArgs>0)
        th=varargin{1};
    end
    in1=0;
    first1=0;
    last1=0;
    for i=1:K
        if in1
            if abs(r(i))<th
                last1=i-1;
                in1=0;
                locs=[locs;last1];            
            end
        else
            if abs(r(i))>=th                            
                first1=i;
                in1=1;
            end
        end        
    end
    if in1
        last1=i;
        locs=[locs;last1];
    end
end


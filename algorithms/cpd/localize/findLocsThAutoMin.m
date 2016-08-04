function locs=findLocsThAutoMin(r,varargin)
    if(size(varargin,1)<1)
        p=0.01;
    else
        p=varargin{1};
    end
    r2=sort(r,'ascend');    
    mx=p*r2(end);
    i=0; s=0;
    while(i<length(r2) && s<mx);    i=i+1;   s=s+r2(i);   end
    if(i==0); th=0; else th=r2(i); end;
    
    K=numel(r);
    locs=[];      
    in1=0;
    first1=0;
    last1=0;
    for i=1:K
        if in1
            if abs(r(i))<th
                last1=i-1;
                in1=0;
                locs=[locs;first1];            
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
        locs=[locs;first1];
    end
end

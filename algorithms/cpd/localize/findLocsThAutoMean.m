function locs=findLocsThAutoMean(r,varargin)
    delta=0.1;
    if(size(varargin,1)<1)
        p=0.01;
    else
        p=varargin{1};
    end
    r2=sort(r,'descend');    
    th=max(eps,r2(min(numel(r2),max(1,round(numel(r2)*p*r2(1)*(1+delta))))));
    %r2=sort(r,'ascend');    
    %mx=p*r2(end);
    %i=0; s=0;
    %while(i<length(r2) && s<mx);    i=i+1;   s=s+r2(i);   end    
    %if(i==0); th=0; else th=r2(i); end;
    
    r=r./max(r);
    th=1e-10;
    
    K=numel(r);
    locs=[];      
    in1=0;
    first1=0;
    last1=0;
    lastval=th;
    for i=1:K
        if in1
            if abs(r(i))<firstval
                last1=i-1;
                in1=0;
                locs=[locs;min(K,max(1,round(0.5*(first1+last1))))];            
            end
        else
            if abs(r(i))>=lastval
                first1=i;
                in1=1;
                firstval=r(i);
            end
        end        
    end
    if in1
        last1=i;
        locs=[locs;min(K,max(1,round(0.5*(first1+last1))))];
    end
end

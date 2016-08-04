function m=matcc(tp,fp,tn,fn)
    d=(tp+fp).*(tp+fn).*(tn+fp).*(tn+fn);    
    m=(tp.*tn-fp.*fn)./(sqrt(d));
    m(d==0)=0;    
end
figure; 
aa=[
    -0.7,0.5,0.9;
    -0.2,0.5,0.9;
    -0.7,0.0,0.09;
    -0.0,0.5,0.04;
    ];
b=[1,2,3,2,1];
for i=1:4; 
    subplot(2,2,i); 
    a=aa(i,:)'; 
    x=generateARMA(a,b,100,1); plot(x); 
    title(sprintf('%3.2f, ',a));  
end;
figure; 
d=[0,1,2,3];
aa=[
    -0.7,0.5,0.9;
    -0.7,0.5,0.9;
    -0.7,0.5,0.9;
    -0.7,0.5,0.9;
    ];
b=[1,2,3,2,1];
for i=1:4; 
    subplot(2,2,i); 
    a=aa(i,:)'; 
    x=generateARIMA(a,b,d(i),100,1); plot(x); 
    title(sprintf('a=[%3.2f,%3.2f,%3.2f],d=%d ',a(1),a(2),a(3),d(i)));  
end;
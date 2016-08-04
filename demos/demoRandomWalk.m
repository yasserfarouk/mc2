figure; 
for i=1:9; 
    subplot(3,3,i); 
    p=0.1*(i-1)./0.8; 
    plot(generateRandomWalk(100,3,1,p)); 
    title(sprintf('p=%f',p));  
end;
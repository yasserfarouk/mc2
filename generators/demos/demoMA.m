figure; 
for i=1:8; 
    subplot(4,2,i); 
    a=[1,repmat(2*ones(1,3),1,i)]; 
    x=generateMA(a,1000,1); plot(x); 
    title(sprintf('$m=%d,\\left\\langle X\\right\\rangle=%5.3f,Var\\left(X\\right)=%5.3f$',numel(a),mean(x),var(x)),'interpreter','latex');  
end;
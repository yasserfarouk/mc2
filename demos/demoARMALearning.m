rng default;
A=[-0.7,0.5,0.9];
B=[1,2,3,2,1];


[x,noise]=generateARMA(A,B,1024);
[a,b,sigma2,x_0_m,a0,mu_epsilon,eps]=learnARMALS(x,3,5);
rng default; 
y=generateARMA(a,b,1024); 


estimateAvailable=true;
try
    ToEstMdl = arima(3,0,5); 
    ToEstMdl.Constant = 0;
    EstMdl = estimate(ToEstMdl,x);
    rng default;
    z=generateARMA(cell2mat(EstMdl.AR),cell2mat(EstMdl.MA),1024);
catch
    estimateAvailable=false;
end

figure;
plot(x); hold on; plot(y); legend({'Original time-series','Predicted time-series'});

if estimateAvailable
figure;
subplot(2,1,1); hold on; plot([x,y,z]); 
legend({'Original time-series','Two-stages Regression','Maximum Likelihood'});

subplot(2,1,2); plot([y-x,z-x]);
ylabel('Error in prediction')
legend({'Two-stages Regression','Maximum Likelihood'});
end

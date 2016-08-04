T=100;
ahat=0.05;
bhat=5;
sigma2=1;
xhat=generateARMA([-0.7,0.5,.9],[1,2,3,2,1],T,1); 
R=sigma2.*randn(T,1);

% A linear trend
x=xhat+ahat.*(1:T)'+bhat+R;
[y,a,b]=detrendLinear(x);
figure; 
subplot(2,2,1); 
plot(x,'LineWidth',2); 
hold on;
plot(ahat.*(1:T)'+bhat,'g','LineWidth',2);
plot(a.*(1:T)'+b,'r','LineWidth',2);
title('$X=ARMA(3,5)+T_0+R, ra_t \sim \mathcal{N} (0,1)$', 'interpreter','latex'); 
legend({'Mixed Time-Series','Linear Trend','Best Fit Line'});
subplot(2,2,3); 
plot(y,'LineWidth',2);
hold on; 
plot(xhat,'r','LineWidth',2);
title('$T_0=0.5t+5.0$', 'interpreter','latex');
legend({'Detrended time-series','Groud Truth'});

% a nonlinear trend

T_0=cumsum(0.5*abs(randn(T,1)));
x=xhat+R+T_0;
[y,a,b]=detrendLinear(x);
subplot(2,2,2); 
plot(x,'LineWidth',2); 
hold on;
plot(T_0,'g','LineWidth',2);
plot(a.*(1:T)'+b,'r','LineWidth',2);
title('$X=ARMA(3,5)+T_0+R , ra_t \sim \mathcal{N} (0,0.5)$', 'interpreter','latex'); 
legend({'Mixed Time-Series','Linear Trend','Best Fit Line'});
subplot(2,2,4); 
plot(y,'LineWidth',2);
hold on; 
plot(xhat,'r','LineWidth',2);
title('$T_0=CUMSUM\left(\left\|R\right\|\right), r_t \sim \mathcal{N} (0,3)$','interpreter','latex');
legend({'Detrended time-series','Groud Truth'});
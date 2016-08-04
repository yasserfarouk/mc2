
figure;
Ns=2;
Nt=1;
sigma2=2;
A=[0.9,0.1;0.1,0.9];
mu=[6,-6];
p0=rand(Ns,1);
p0=p0./sum(p0);

sigma=zeros(Nt,Nt,Ns);
letters={'(a)','(b)','(c)','(d)'};
for i=1:Ns
    sigma(:,:,i)=sigma2*eye(Nt);
end

xhmm=generateHMM(A,p0,mu,sigma,100);
[p0L,AL,muL,sigmaL,log_pX,err,iter]=learnHMM(xhmm,Ns);
xL=generateHMM(AL,p0L,muL,sigmaL,100);

% figure;
% plot(xhmm,'LineWidth',2);
% hold on;
% plot(xL,'LineWidth',2);
% legend({'Original','Learned'});

Nseries=5;
xs=cell(Nseries,1);
Ts=[100,120,150,130,80];
for s=1:Nseries
    xs{s}=generateHMM(A,p0,mu,sigma,Ts(s));
end
[p0L,AL,muL,sigmaL]=learnHMMMulti(xs,Ns);

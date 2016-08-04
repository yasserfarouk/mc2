T=3000;
E=4;
L=[10:10:T,T];
useHankelization=true;

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 20)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 20)

set(0,'defaultlinelinewidth',2)


x=generateCoupled('two',[0.4,0.2,3.8,3.5,1,0],T);

rho=zeros(numel(L),2);
for l=1:numel(L)
    fprintf('%d of %d (%d)\n',l,numel(L),L(l));
    [g,y,gDetails]=ccm(x,[],E,L(l),useHankelization);
    rho(l,1)=g(1,2);
    rho(l,2)=g(2,1);
end

figure;
hold on;
plot(L,rho(:,1));
plot(L,rho(:,2));
ylabel('\rho');
xlabel('T');
legend({'corr(X, X|M_y)','corr(Y, Y|M_x)'});
set(gcf,'units','normalized','outerposition',[0 0 1 1]);    
savefig(strcat('rho01large','.png'), 'png', '-rgb', '-c0.1', '-r250');

figure;
hold on;
if useHankelization
    plot(x(:,1)-squeeze(y(1,2,:))); title('Error in predicting x');
else
    plot(x(E:end,1)-squeeze(y(1,2,:))); title('Error in predicting x');
end
set(gcf,'units','normalized','outerposition',[0 0 1 1]);    
savefig(strcat('x01large','.png'), 'png', '-rgb', '-c0.1', '-r250');

figure;
hold on;
if useHankelization
    plot(x(:,2)-squeeze(y(2,1,:))); title ('Error in prdicting y');
else
    plot(x(E:end,2)-squeeze(y(2,1,:))); title ('Error in predicting y');
end
set(gcf,'units','normalized','outerposition',[0 0 1 1]);    
savefig(strcat('y01large','.png'), 'png', '-rgb', '-c0.1', '-r250');


figure;
hold on;
plot(x(1:100,1),'b');
plot(x(1:100,2),'r--');
legend({'x','y'});
set(gcf,'units','normalized','outerposition',[0 0 1 1]);    
mysavefig(strcat('xy01large','.png'), 'png', '-rgb', '-c0.1', '-r250');

%[M,G]=detectCCMC(x',0.8,4);
%G
data=dlmread('./Data/simInfoIn05.txt');     % this file can be generated using a script in imitationNav project
motors=[zeros(20,2);data(:,2:3)]; state=[zeros(20,3);data(:,4:6)];

subplot(2,1,1);
pmotors=tspca(motors','norm',1,'plot',1);
title('motor command after normalization');

subplot(2,1,2);
pstate=tspca(state','norm',1,'plot',1);
title('state after normalization');

figure;
subplot(3,1,1);
plot(pmotors,'r'); hold on;
plot(pstate,'b');hold off;
title('motor commands and state projection');

[rmotors,lmotors]=rsst(pmotors,5,6,-1,'th',0.01); %,'dfunname','distWeightedEigDist');
[rstate,lstate]=rsst(pstate,15,25,-1,'th',0.01); %,'dfunname','distWeightedEigDist');

subplot(3,1,2);
plot(rmotors,'r'); hold on;
plot(rstate,'b'); hold off;
title('RSST output (change scores)');

locm=zeros(1,size(motors,1)); locm(lmotors)=1;
locs=zeros(1,size(state,1)); locs(lstate)=1;

subplot(3,1,3);
plot(locm,'r'); hold on;
plot(locs,'b'); hold off;
title('Localization output (change locations)');


n=numel(pmotors);
pmotors(2:n)=abs(pmotors(2:n)-pmotors(1:n-1));
pstate(2:n)=abs(pstate(2:n)-pstate(1:n-1));
pmotors(1)=0;
pstate(1)=0;

figure;
subplot(3,1,1);
plot(pmotors,'r'); hold on;
plot(pstate,'b');hold off;
title('motor commands and state projection');

[rmotors,lmotors]=rsst(pmotors,15,25,-1,'th',0.01); %,'dfunname','distWeightedEigDist');
[rstate,lstate]=rsst(pstate,15,25,-1,'th',0.01); %,'dfunname','distWeightedEigDist');

subplot(3,1,2);
plot(rmotors,'r'); hold on;
plot(rstate,'b'); hold off;
title('RSST output (change scores) OF DIFF');

locm=zeros(1,size(motors,1)); locm(lmotors)=1;
locs=zeros(1,size(state,1)); locs(lstate)=1;

subplot(3,1,3);
plot(locm,'r'); hold on;
plot(locs,'b'); hold off;
title('Localization output (change locations)  OF DIFF');



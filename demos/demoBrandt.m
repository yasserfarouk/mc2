N=200;
L=100;
T=3000;
range=[-5,5];    
scale2range=1;
noiseLevel=0.00;
gNoiseSigmaInitial=1.0;       % initial noise
pOutlier=0.0;
minChangeDistance=150;
pchange=0.01;

ps=[5,10,30];
qs=[0];

[xbase,locsTbase,tbase,xTruebase]=produceSingle(T,'pchange',pchange,...
    'minChangeDistance',minChangeDistance,'noiseLevel',noiseLevel,...
    'poutlier',pOutlier,'noiseSigma',gNoiseSigmaInitial,'range',range...
    ,'scale2range',scale2range);

c=zeros(T,numel(ps)*numel(qs));
labels=cell(numel(ps)*numel(qs),1);
for i=1:numel(ps)
    for j=1:numel(qs)
        fprintf('%d,%d\n',i,j);
        c(:,(i-1)*numel(qs)+j)=brandt(xbase,N,L,ps(i),ps(i)*qs(j));
        labels{(i-1)*numel(qs)+j}=sprintf('p=%d,q=%d',ps(i),ps(i)*qs(j));
    end
end

figure;
subplot(2,1,1);
plot(xbase);
subplot(2,1,2);
plot(c);
legend(labels,'Location','northoutside','Orientation','horizontal');
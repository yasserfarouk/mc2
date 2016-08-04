groundTruth=sin([1:500]*pi/50)+2*sin([1:500]*pi/5)+[-249:250]./50;
noise=randn(1,500);
x=groundTruth+noise;

imf=emd(x);
n=size(imf,1);
nRow=ceil(n/2)+1;
nCol=2;
figure;
subplot(nRow,nCol,1); plot(groundTruth,'LineWidth',2); title('Ground Truth', 'FontSize',20);
subplot(nRow,nCol,2); plot(x,'LineWidth',2); title('Input after adding Gaussian noise', 'FontSize',20);
for i=1:n
    subplot(nRow,nCol,i+2);
    plot(imf(i,:),'LineWidth',2);
    title(sprintf('IMF #%d',i),'FontSize',20);
end




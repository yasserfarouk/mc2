close all;
y=zeros(2,400);
y(1,:)=[1.5*sin(.2.*(1:200)),1.5*sin(.3.*(1:200))];
plot(y(1,:));
s=1;
p=1;p2=1;
f1=figure;
f2=figure;
y=y+randn(2,400);
plot(y(s,:));
for i=0:1
    for j=0:1
        for k=0:1
            [r,c]=sst3(y(s,:),40,41,-1,-1,-1,-1,-1,'nobalance',i,'normalize',j,'autothreshold',k,'lfuncalconce',-1);            
            %[r,c,rloc,cloc]=rsst(y(s,:),40,41,-1,'nobalance',i,'normalize',j,'autothreshold',k,'lfuncalconce',-1);            
            figure(f2);
            subplot(4,2,p); 
            plot(c,'r');
            xlabel(sprintf('CUMSUM,nonbalance %d,normalize,%d,autothreshold,%d',i,j,k));                        
            p=p+1;
        end
        figure(f1);
        subplot(2,2,p2); 
        plot(r,'b');
        xlabel(sprintf('CP,nonbalance %d,normalize,%d,autothreshold,%d',i,j,k));
        p2=p2+1;
    end
end

% noiseLevel=0.01;
% gNoiseSigma=0.1;
% w=5;
% n=10;
% 
% [x,motifsT,locsT,occursT,xTrue]=generatePatterns(3,[3,3,4],[50,200],1,[-noiseLevel/2,noiseLevel/2],100,noiseLevel,gNoiseSigma);
% t=zeros(1,numel(x));
% for i=1:length(locsT)
%     r=locsT{i};
%     for j=1:size(r,1)
%         t(r(j,1))=1;
%         t(r(j,2))=1;
%     end
% end
% 
% 
% ch=[];
% times=[];
% i=1;
% chNames=cell(1,1);
% 
% tmp=tic;
% y=t;
% tmp=toc(tmp); times=[times tmp];
% ch=[ch;y];
% ch(size(ch,1),:)=ch(size(ch,1),:)./max(ch(size(ch,1),:));
% chNames{i,1}='True'; i=i+1;
% 
% tmp=tic;
% y=sst3(x,w,n,3,0,n);
% tmp=toc(tmp); times=[times tmp];
% ch=[ch;y];
% ch(size(ch,1),:)=ch(size(ch,1),:)./max(ch(size(ch,1),:));
% chNames{i,1}='0,n'; i=i+1;
% 
% % tmp=tic;
% % y=sst2(x,w,n,3,0,n);
% % tmp=toc(tmp); times=[times tmp];
% % ch=[ch;y];
% % ch(size(ch,1),:)=ch(size(ch,1),:)./max(ch(size(ch,1),:));
% % chNames{i,1}='SST (HH'')'; i=i+1;
% 
% tmp=tic;
% y=sst3(x,w,n,3);
% tmp=toc(tmp); times=[times tmp];
% ch=[ch;y];
% ch(size(ch,1),:)=ch(size(ch,1),:)./max(ch(size(ch,1),:));
% chNames{i,1}='def'; i=i+1;
% 
% tmp=tic;
% y=sst3(x,w,n,3,0,n,0,n);
% tmp=toc(tmp); times=[times tmp];
% ch=[ch;y]; 
% ch(size(ch,1),:)=ch(size(ch,1),:)./max(ch(size(ch,1),:));
% chNames{i,1}='0,n,0,n'; i=i+1;
% 
% 
% times=times./(numel(x)/1000);       % find time per point in ms
% 
% n=size(ch,1);
% subplot(n+2,1,1); plot(x); hold on; plot(t,'r'); hold off;
% subplot(n+2,1,2); plot(xTrue); hold on; plot(t,'r'); hold off;
% for i=1:n
% subplot(n+2,1,i+2); plot(ch(i,:)); ylabel(chNames{i,1}); hold on; plot(t,'r'); hold off;
% end
% for i=1:n
%     ch(i,:)=ch(i,:)./max(ch(i,:));
% end
% 
% cstring='rgbcmyk';
% 
% 
% fp=zeros(11,n);
% tp=zeros(11,n);
% tn=zeros(11,n);
% fn=zeros(11,n);
% mc=zeros(11,n);
% the=0.1;
% for i=1:11
%     th=(i-1)*the;
%     [tp(i,:),fp(i,:),tn(i,:),fn(i,:)]=cpquality(ch',t,w,th);
%     mc(i,:)=matcc(tp(i,:),fp(i,:),tn(i,:),fn(i,:));
% end
% 
% pos=sum(t);
% acc=(tp+tn)./numel(t);
% sp=tp./(pos);
% prec=tp./(tp+fp);
% recall = tp./(tp+fn);
% F=2.*(prec.*recall)./(prec+recall);
% 
% figure;
% hold on;
% for i=1:n
% plot(prec(:,i),recall(:,i),cstring(mod(i,7)+1)); 
% end
% title('recall vs. prec');
% legend(chNames); 
% 
% 
% figure;
% hold on;
% for i=1:n
% plot(mc(:,i),cstring(mod(i,7)+1)); 
% end
% title('Matthews correlation coefficient');
% legend(chNames); 
% 
% 
% figure;
% hold on;
% for i=1:n
% plot((0:10).*the,F(:,i),cstring(mod(i,7)+1)); 
% end
% title('F-measure vs. threshold');
% legend(chNames); 
% 
% 
% figure;
% hold on;
% for i=1:n
% plot(acc(:,i),sp(:,i),cstring(mod(i,7)+1)); 
% end
% title('specificity vs. accuracy');
% legend(chNames); 
% 
% figure;
% hold on;
% for i=1:n
% plot(fp(:,i),tp(:,i),cstring(mod(i,7)+1)); 
% end
% title('ROC (tp vs. fp)');
% legend(chNames); 
% % 
% nr=round(numel(t)/100);
% q=zeros(nr,n);
% for i=1:nr    
%     [q(i,:)]=cpsp(ch',t,i);
% end
% 
% figure;
% hold on;
% for i=1:n
% plot(q(:,i),cstring(mod(i,7)+1)); 
% end
% title('My Change Point Specificity');
% legend(chNames,'Location','NorthWest'); 
% hold off;
% 
% 
% h=figure;
% axes1 = axes('Parent',h,'XTickLabel',chNames);
% box(axes1,'on');
% hold(axes1,'all');
% plot(times,'-or');
% 

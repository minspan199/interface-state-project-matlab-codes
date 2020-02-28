clear all;

ksA = 0;t1 = 0.5;t2 = 1;k = 1;n = 15;
H = diag([repmat([t2 t1],1,n) t2],1) + diag([repmat([t2 t1],1,n) t2],-1);
[V,D] = eig(H);
figure;plot(diag(D),'bo');hold on;
set(gca,'Fontname','Times New Roman');set(gcf, 'Position', [00, 00, 350, 200])
set(gca,'Fontsize',14);xlim([0 30]);

figure;bar(abs(V(:,16))/max(abs(V(:,16))),'b');hold on;
set(gca,'Fontname','Times New Roman');set(gcf, 'Position', [00, 00, 400, 200])
set(gca,'Fontsize',14);xlim([0 30])

H = diag([repmat([t2 t1],1,7) repmat([t1 t2],1,7)],1) + diag([repmat([t2 t1],1,7) repmat([t1 t2],1,7)],-1);
[V,D] = eig(H);
figure;plot(diag(D),'bo');hold on;
set(gca,'Fontname','Times New Roman');set(gcf, 'Position', [00, 00, 350, 200])
set(gca,'Fontsize',14);xlim([0 30]);
figure;bar(abs(V(:,15))/max(abs(V(:,15))),'b');hold on;
xticks([0 5 10 15 20 25 30]);xticklabels({'-15','-10','-5','0','5','10','15'});
set(gca,'Fontname','Times New Roman');set(gcf, 'Position', [00, 00, 300, 200])
set(gca,'Fontsize',14);xlim([0 30])
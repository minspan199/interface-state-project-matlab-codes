clear all;
close all;
beta = 1;
kapa = 1;
kama = 1;
figure;
hold on;
for k = -pi/2:0.01:pi/2
    H=[beta kapa 0 -kama*exp(-1i*k);
        kapa beta kama*exp(1i*k) 0;
        0 kama*exp(-1i*k) -beta -kapa;
        -kama*exp(1i*k) 0 -kapa -beta];
    lam = eig(H);
    [~,ind]=max((real(lam)));
    plot(k, lam(ind), 'b.')
    [~,ind2]=min((real(lam)));
    plot(k, lam(ind2), 'b.')
    lam_p=find(lam>0);
    ind = min(real(lam_p));
    plot(k, lam(ind), 'r.')
    lam_m=find(lam<0);
    ind = max(real(lam_m));
    plot(k, lam(ind), 'r.')
end

set(gca,'Fontname','Times New Roman');
set(gcf, 'Position', [00, 00, 175, 150])
set(gca,'Fontsize',14);
axis([-pi/2 pi/2 -3 3]);
box on;
set(gca,'xtick',[-pi/2,0,pi/2]);
set(gca,'xticklabel',{'-\pi/2','0','\pi/2'});


beta_ = beta^2 + kama^2 + kapa^2;
kama_ = kama*kapa;
kapa_ = 2*beta*kapa;
figure;
hold on;
for k = -pi:0.01:pi
    H = [beta_ + 2*kama_*cos(k) kapa_;
        kapa_ beta_ - 2*kama_*cos(k)];
    lam1 = eig(H);
    
    plot(k, lam1, 'black.')

end

axis([-pi pi 0 6]);
box on;
set(gca,'xtick',[-pi,-pi/2,0,pi/2,pi]);
set(gca,'xticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});
set(gca,'Fontname','Times New Roman');
set(gcf, 'Position', [00, 00, 175, 150])
set(gca,'Fontsize',14);



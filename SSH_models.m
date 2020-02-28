clear all;

ksA = 0;
t1 = 0.5;t2 = 1;
k = 1;


r = 0;
for kb = -pi:0.05:pi
    H = [1i*r t1 + t2*exp(-1i*kb);t1 + t2*exp(1i*kb) -1i*r];
    tem{k} = eig(H);
    ksp(k) = tem{k}(1);
    ksm(k) = tem{k}(2);
    k = k+1;
end

figure

plot(-pi:0.05:pi,real(ksp),'b.');hold on;
% plot(-pi:0.05:pi,imag(ksp),'r.');
plot(-pi:0.05:pi,real(ksm),'r.');
% plot(-pi:0.05:pi,imag(ksm),'r.');
set(gca,'xtick',[-pi,-pi/2,0,pi/2,pi]);
set(gca,'xticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});
set(gca,'Fontname','Times New Roman')
set(gcf, 'Position', [00, 00, 250, 200])
set(gca,'Fontsize',14)
% legend('Re[\epsilon]','Im[\epsilon]')
xlim([-pi pi]);ylim([-2 2])


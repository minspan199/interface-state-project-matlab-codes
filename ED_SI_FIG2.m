clear all;

ksA = 0;
ta = 0.031;
tb = 0.016;
k = 1;
kk = 1;
r=0;
% r = 0:0.01:1

figure
for r = 0.05:0.05:0.05
    k = 1;
    for kb = -pi:0.05:pi
        H = [1i*r ta + tb*exp(-1i*kb);ta+tb*exp(1i*kb) -1i*r];
        tem{k} = eig(H);
        ksp(k) = tem{k}(1);
        ksm(k) = tem{k}(2);
        k = k+1;
    end
    kb = -pi:0.05:pi;
    subplot(2,1,1)
    plot(kb,real(ksp),'b.');
    set(gca,'xtick',[-pi,-pi/2,0,pi/2,pi]);
    set(gca,'xticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});
    h = get(gca, 'xlabel');
    set(h, 'FontName', 'Arial')
    h = get(gca, 'ylabel');
    set(h, 'FontName', 'Arial')
    tex.Fontname = 'Arial';
    axis([-pi pi -0.05 0.05])
    set(gca,'Fontname','Arial')
    set(gca, 'FontSize', 12)
%     set(gca,'Fontsize',12)
    hold on
    plot(kb,real(ksm),'b.');
%     xlabel('ka','FontSize',14,'FontName','Arial');
%     ylabel('Re[\epsilon]','FontSize',14,'FontName','Arial');
    hold on
    
    %%%%%%%%%%%%%%%%%%
    subplot(2,1,2)
    plot(kb,imag(ksp),'b.');
    hold on
    plot(kb,imag(ksm),'b.');
    set(gca,'xtick',[-pi,-pi/2,0,pi/2,pi]);
    set(gca,'xticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});
    set(gca, 'FontSize', 12)
    if kk == 1
        axis([-pi pi -0.05 0.05])
    end
    if kk == 2
        axis([-pi pi -0.05 0.05])
    end
    hold on
    set(gca,'xtick',[-pi,-pi/2,0,pi/2,pi]);
    set(gca,'xticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});
    set(gca,'Fontname','Arial')
%     set(gca,'Fontsize',12)
    hold on
%     xlabel('ka','FontSize',14,'FontName','Arial');
%     ylabel('Im[\epsilon]','FontSize',14,'FontName','Arial');
    kk = kk+1;
end
set(gcf, 'Position', [00, 00, 300, 300])
% figure
% plot(imag(ksm))
% plot(imag(ksp))
% hold on
% plot(imag(ksm))
% 
% k=1;   
% 
% global r
% for R = 0:0.0005:0.1
% r=R;
%     xx(k) = R;
%    pk = abs(ta+tb*exp(1i*ka));
%    yy(k) = quad(@myfun0,-2*pi,2*pi);
%    k=k+1;
% end
%     
% figure
% plot(xx,real(yy));
% hold on
% plot(xx,imag(yy));
% 
% function y = myfun0(ka) 
% global r
% r
% ta = 0.031;
% tb = 0.016;
%    y = 1/2*(1+(1i*r./sqrt((abs(ta+tb*exp(1i*ka))).^2-r.^2)));
% end
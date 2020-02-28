clear all
close all
D = -70:1:250; % distances of the waveguide shifts.
Ta1 = [0.0301 0.0283 0.0267 0.0251 0.0236 0.0222 0.0209 0.0197 0.0185 0.0174 0.0164 0.0154 0.0145 0.0136 0.0128 0.0121 0.0114 0.0107 0.01 0.0095 0.0089 0.0084 0.0079 0.0074 0.007 0.0066 0.0062 0.0058 0.0055 0.0052 0.0049 0.0046 0.0043]; % coupling parameters
Tb1 = [0.0052 0.0055 0.0058 0.0062 0.0066 0.0070 0.0074 0.0079 0.0084 0.0089 0.0095 0.0100 0.0107 0.0114 0.0121 0.0128 0.0136 0.0145 0.0154 0.0164 0.0174 0.0185 0.0197 0.0209 0.0222 0.0236 0.0251 0.0267 0.0283 0.0301 0.0320 0.0341 0.0362]; % coupling parameters
for kk = 0:1:31
    Ta(10*kk+1) = Ta1(kk+1);
    d = (Ta1(kk+1)-Ta1(kk+2))/9;
    Ta(10*kk+2:10*kk+10) = Ta1(kk+1)-d:-d:Ta1(kk+2);
    
    Tb(10*kk+1) = Tb1(kk+1);
    d = (Tb1(kk+1)-Tb1(kk+2))/9;
    Tb(10*kk+2:10*kk+10) = Tb1(kk+1)-d:-d:Tb1(kk+2);
    
end  % I tested the coupling at 33 points, here to draw the continuum, the parameters were interpolated to 321 points.
Ta(321) = Ta1(33);Tb(321)=Tb1(33); % coupling parameter at 250 nm of shifting
figure
hold on
for k = 1:1:321
 
m = 25;
n = 25;
DD1 = [1i*0.02 1i*0.08]; % Loss contrast in the left lattice. 
b1 = repmat(DD1,[1 m]); 
DD2 = [0 1i*0.02];  % Loss contrast in the right lattice. 
b2 = repmat(DD2,[1 n]);
 
t1 = [repmat([0.008 0.02 ],[1 m-1]) 0.008 Ta(k) Tb(k) repmat([0.02 0.008],[1 n-1])];
H = diag([b1 b2])+diag(t1,1)+diag(t1,-1); 
[V{k}, A] = eig(H);
lam{k} = diag(A);
[a(k), b(k)] = min(abs(lam{k}));
lan(k) = lam{k}(b(k));
B(k,:) = real(lam{k}((real(lam{k}))>0.0001));
C(k,:) = real(lam{k}((real(lam{k}))<-0.0001));
plot(D(k),B(k,:),'blue.')
hold on
plot(D(k),C(k,:),'blue.')
end
 
hold on
plot(D,real(lan),'LineWidth', 2)
set(gcf, 'Position', [00, 00, 300, 250])
% xlabel('Waveguide shifting D (nm)','FontSize',14,'FontName','Arial');
% ylabel('Effective index','FontSize',14,'FontName','Arial');
axis([-70 250 -0.045 0.045])
box on


figure
plot((-49:1:50),(abs(V{71}(:,b(71)))./max(abs(V{71}(:,b(71))))));

hold on
% xlabel('Waveguide Site Number','FontSize',14,'FontName','Arial');
% ylabel('Intensity','FontSize',14,'FontName','Arial');
set(gcf, 'Position', [00, 00, 300, 250])
set(gcf, 'Position', [00, 00, 320, 160])
set(gca,'FontSize',12)
xticks([-50 -25 0 25 50])
xticklabels({'-50','-25','0','25','50'})
yticks([0 0.2 0.4 0.6 0.8 1])
yticklabels({'0' '0.2' '0.4' '0.6' '0.8' '1'})

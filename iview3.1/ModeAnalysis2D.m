close all
clc
clear all
t = 1;
% N = 8;
figure
r = 2;
% mov = VideoWriter('Movie0.avi');
% mov.FrameRate = 1;
% open(mov)

% D = -0:-0.005:-3;
k = 1;
% D = -1;
D = diag([0 1i*r 0 1i*r])+diag(t*ones(1,3),1)+diag(t*ones(1,3),-1);
E = fliplr(diag(1*t*ones(1,4)));
F = zeros(4);
H = [D E F F;E D E F;F E D E;F F E D];
T = 1;


while 1
 

p = zeros(1,16);
p(1) = -1i*T;
H = H + diag(p);

[A,V] = eig(H);
Re = diag(real(V)/t);
Im = diag(imag(V)/t);

plot(Re(14),Im(14),'b*');
plot(Re(13),Im(13),'r*');
% axis([-2 2 -3 4])
hold on
% text(80,220,strcat('D = ', num2str(D(k))),'FontSize',16)
yPos = 0;
hold on
plot(get(gca,'xlim'), [yPos yPos],'b--'); % Adapts to x limits of current axes

xlabel('Real')
ylabel('Imaginary')
% TEM = getframe;
% writeVideo(mov,TEM.cdata);
lam = diag(V);
if find(imag(lam)<-1e-6)
    break;
else
    T = T+0.002*t;
end
k = k+1;
end
% close(mov);

b = find(imag(diag(V))<-1e-6);
% b=4;
M = abs(A(:,b)).*abs(A(:,b));
M = M./max(M);
N = angle(A(:,b));

figure;
bar(M)
ylim([0 1]);
figure;
bar(N)
ylim([-pi pi])
set(gca,'ytick',[-pi,-pi/2,0,pi/2,pi]);
set(gca,'yticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});

RT = [M(1:4)';fliplr(M(5:8)');M(9:12)';fliplr(M(13:16)')];
RA = [N(1:4)';fliplr(N(5:8)');N(9:12)';fliplr(N(13:16)')];

figure
imagesc(RT);
colormap jet
colorbar
figure
imagesc(RA);
colormap([0 0 0;1 0 0;0 1 0;0 0 1;0 0 0]);
colorbar% X(:,2) = A(:,4);
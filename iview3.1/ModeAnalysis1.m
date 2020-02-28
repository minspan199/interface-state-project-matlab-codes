close all
clear all
t = 1;
N = 8;
figure
r = -2;
mov = VideoWriter('Movie0.avi');
mov.FrameRate = 1;
open(mov)

D = -0:0.005:5;
k = 1;
% D = -1;
while 1
T = diag([0 1i*r 0 1i*r 0 1i*r 1i*D(k) 1i*r]);
% T1 = diag(t*ones(1,N-1),1);
% T2 = diag(t*ones(1,N-1),-1);

T1 = diag([t 0.75*t t 0.75*t t 0.75*t t],1);
T2 = diag([t 0.75*t t 0.75*t t 0.75*t t],-1);
H = T + T1 + T2;
[A,V] = eig(H);
Re = diag(real(V)/t);
Im = diag(imag(V)/t);
% hold off
plot(Re(1), Im(1), 'c*')    %   Plot real and imaginary parts
hold on
plot(Re(2), Im(2), 'c.','LineWidth',2)    %   Plot real and imaginary parts
plot(Re(3), Im(3), 'b*')    %   Plot real and imaginary parts
plot(Re(4), Im(4), 'b.','LineWidth',2)    %   Plot real and imaginary parts
plot(Re(5), Im(5), 'r*')    %   Plot real and imaginary parts
plot(Re(6), Im(6), 'r.','LineWidth',2)    %   Plot real and imaginary parts
plot(Re(7), Im(7), 'g*')    %   Plot real and imaginary parts
plot(Re(8), Im(8), 'g.','LineWidth',2)    %   Plot real and imaginary parts
% plot(Re(9), Im(9), 'c.')    %   Plot real and imaginary parts
% axis([-2 2 -3 4])
hold on
text(80,220,strcat('D = ', num2str(D(k))),'FontSize',16)
yPos = 0;
hold on
plot(get(gca,'xlim'), [yPos yPos],'b--'); % Adapts to x limits of current axes

xlabel('Real')
ylabel('Imaginary')
% TEM = getframe;
% writeVideo(mov,TEM.cdata);
if find(Im>0)
    break;
end
k = k+1;
end
close(mov);

b = find(Im>0);
X(:,1) = abs(A(:,b)).*abs(A(:,b));
X(:,1) = X(:,1)./max(X(:,1));
X(:,2) = angle(A(:,b));

figure;
bar(X(:,1))
ylim([0 1]);
figure;
bar(X(:,2))
ylim([-pi pi])
set(gca,'ytick',[-pi,-pi/2,0,pi/2,pi]);
set(gca,'yticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});

% X(:,2) = A(:,4);
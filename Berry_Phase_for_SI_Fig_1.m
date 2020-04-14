clear all
% close all
ta = 2.5;
tb = 1;

for kk = 1:1:5000
    kk
    k = -pi:0.00001:pi;
    r = 0.002*(kk-1);
    rouk = ta + tb*exp(-1i*k);
    lam = sqrt((ta + tb*cos(k)).^2 + (tb*sin(k)).^2 - r^2);
    rou = sqrt((ta + tb*cos(k)).^2 + (tb*sin(k)).^2);
    phi = angle(rouk);
    rk = atan(rou/1i/r);
    SIN = rou./lam;
    COS = 1i*r./lam;
    COS2 = -2*(r.^2)/r/r./lam;
    SIN2 = SIN.^2;
    Re1(kk) = inte(phi,COS)/2;
    Re2(kk) = inte(phi,-COS)/2;
end
R = 0.002*((1:1:5000)-1);

figure
plot(R/ta,real(Re1),'b.')
hold on
plot(R/ta,real(Re2),'r.')
axis([0 4 -2 2])
set(gca,'Fontname','Times New Roman');
set(gcf, 'Position', [00, 00, 450, 300])
set(gca,'Fontsize',14);

figure
plot(R/ta,imag(Re1),'b.')
hold on
plot(R/ta,imag(Re2),'r.')
axis([0 4 -2 2])
set(gca,'Fontname','Times New Roman');
set(gcf, 'Position', [00, 00, 450, 300])
set(gca,'Fontsize',14);

figure
subplot(2,2,1)
plot(R/ta,real(Re1),'b.')
hold on
plot(R/ta,real(Re2),'r.')
axis([0 4 -2 2])

subplot(2,2,2)
plot(R/ta,imag(Re1),'b.')
hold on
plot(R/ta,imag(Re2),'r.')
axis([0 4 -2 2])

subplot(2,2,3)
plot(R/ta,real(Re1+Re2),'b.')
hold on
plot(R/ta,real(Re1+Re2),'r.')
axis([0 4 -2 2])

subplot(2,2,4)
plot(R/ta,imag(Re1+Re2),'b.')
hold on
plot(R/ta,imag(Re1+Re2),'r.')
axis([0 4 -2 2])

function I = inte(Int,Val)
N = length(Val);
M = 0;
for k = 1:1:N-1
    M = M + ((Int(k+1)-Int(k)))*(Val(k)+Val(k+1))/2;
end
I = M;
end
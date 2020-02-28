clc
clear all
close all
figure
hold on

tb = 0.008/0.008;
ta = 2.5;
delta = 2.5;sigma = 2.5;
m = 10; n = 10; %50 dimmers in the left lattice, and 50 dimmers in the right lattice.
k =1;

% for delta = 0.1:0.01:0.11
    
for N = 127:1:400
%     N = 1000
    Delta(N) = 0.025*N*tb;
    DD1 = [-1i*delta -1i*Delta(N)]-1i*0.0; % Loss contrast in the left lattice.
    b1 = repmat(DD1,[1 m]); 
    DD2 = [1i*0.0 -1i*sigma]-1i*0.0;  % Loss contrast in the right lattice. 
    b2 = repmat(DD2,[1 n]);
t1 = [repmat([tb ta],[1 m+n-1]) tb]; % coupling parameters

H = diag([b1 b2])+diag(t1,1)+diag(t1,-1);  % Assembling of Hamiltonian matrix
[V, A] = eig(H);
lam = diag(A);
[a, b] = min(abs(lam));
[e, f] = max(abs(V(2*m+1,:)));
[t, h] = max(abs(V(2*m,:)));
% min(imag(lam));
% max(abs(V(2*m,:)));
lan(N) = lam(f); % The energy of the interface state
g = find(abs(real(lam+lam(f)))<1e-10);


if(abs(real(lan(N)))>0.000001)
    laf(N) = lam(g);
else
    tem = size(g);
    if(tem(1)>4)
%         [e, f] = max(abs(V(2*m,:)));
        laf(N) = 1i*Inf;
    else
        Ind = find((abs(V(2*m+1,g))<=abs(V(2*m,g))).*(abs(V(2*m,g))>0.1));
        laf(N) = lam(g(Ind));
%         TEM =N;
    end
end
    
% B = (lam((real(lam))>0.000)); C = (lam((real(lam))<-0.000));% finding bulk states and remove the flat band in PT breaking lattice
% plot(B,'blue.')
% hold on
% plot(C,'blue.')
% hold on
% plot(lan(N),'red.')
% hold on
% plot(laf(N),'blue.')
% hold on
% if imag(lan(N))>-1.25
%     part(k) = V(2*m+1,f)/sum(abs(V(:,f)));
%     Th(k) = lan(N)
%     k = k+1
%     break
% end

end
% end 
% laf(N)  = 0;
% laf(500:634) = Inf+1i*Inf;
% lan(500:634) = Inf+1i*Inf;
figure;plot(Delta,real(lan),'b.');hold on;plot(Delta,real(laf),'b.');
% ylabel('Real spectra','FontSize',14,'FontName','Arial');
% xlabel('Delta/tB')
set(gcf, 'Position', [00, 00, 600, 200]);axis([2.5 10 -2 2]);set(gca,'FontSize',14);


figure

% Delta1 = 0.005*(1:1:TEM)*tb;
plot(Delta,imag(lan),'b.');hold on;plot(Delta,imag(laf),'b.');
% ylabel('Imaginary spectra','FontSize',14,'FontName','Arial');
% xlabel('Delta/tB')
set(gcf, 'Position', [00, 00, 600, 200]);axis([2.5 10 -6 0]);set(gca,'FontSize',14);
figure
[e, f] = max(abs(V(2*m+1,:)));bar(abs(V(:,f))/max(abs(V(:,f))),'Linewidth',1); %plot intensity distribution
xlabel('Real spectra','FontSize',14,'FontName','Arial');ylabel('Imaginary spectra','FontSize',14,'FontName','Arial');
set(gcf, 'Position', [00, 00, 350, 300]);box on

figure;[e, b] = max(abs(V(2*m,:)));
plot(abs(V(:,b))/max(abs(V(:,b))),'Linewidth',1); %plot intensity distribution
hold on;xlabel('Waveguide Site Number','FontSize',14,'FontName','Arial');
ylabel('Intensity','FontSize',14,'FontName','Arial');set(gcf, 'Position', [00, 00, 350, 300])

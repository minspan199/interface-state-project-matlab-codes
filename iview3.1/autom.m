function [Rxx]=autom(x)
% [Rxx]=autom(x)
% This function Estimates the autocorrelation of the sequence of
% random variables given in x as: Rxx(1), Rxx(2),…,Rxx(N),
% where N is
% Number of samples in x.
N=length(x);
Rxx=zeros(1,N);
for m=1: N+1
for n=1: N-m+1
Rxx(m)=Rxx(m)+x(n)*x(n+m-1);
end
end
end
function [x,s]=semispanningsdrain(q,L,n,kD,c);
% drain in een semi spanningspakket
% TO 001006
%

dL=L/n;
Lambda=sqrt(kD*c);
x =[-L/2:dL:L/2];
xm=0.5*(x(1:end-1)+x(2:end));

s=zeros(size(x));
for i=1:length(xm)
   rl=abs((x-xm(i))/Lambda);
   s=s+besselk(0,rl);
end
s=s*q*L/(n*2*pi*kD);

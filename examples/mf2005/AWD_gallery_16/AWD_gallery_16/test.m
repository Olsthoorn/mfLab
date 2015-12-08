cL=10;
c0=1;L=100;
x=1:L;
lambda=L*(-1:0.1:1);
clr=repmat('brgkmcy',[1,10]);
 figure; hold on
leg{length(lambda)}='leg';
for i=1:length(lambda)
    c1=(c0-cL*exp(-L/lambda(i)))/(1-exp(-L/lambda(i)));
    c=c1+(cL-c1)*exp(-(L-x)/lambda(i));
    plot(x,c,clr(i)); leg{i}=sprintf('lambda = %g',lambda(i));
end
set(gca,'ylim',[0 cL]);
legend(leg,2);
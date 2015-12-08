% test Boulton 1963

t=logspace(-3,3,60);
r=logspace(-2,3,50);
alpha=0.5;
Sa=1e-5;
Sy=0.1;
he=20; z=25;
t=1;
r=100;
kv=1;
B=15;
T=1000;

Y=boulton2(t,r,T,alpha,Sa,Sy);

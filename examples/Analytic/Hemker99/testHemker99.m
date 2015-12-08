Q = [  1  1  1  1]';
D = [10 10 10 10]';
kv=[  1  1   1 1]';
kh=[ 10 10 10 10]';
Ss=[  1  1  1  1]'*1e-5;
r=10;
rw=0.25;
rc=1;

DD = hemker99(Q,rw,rc,r,D,t,Ss,kh,kv);

figure;
plot(t',squeeze(DD));

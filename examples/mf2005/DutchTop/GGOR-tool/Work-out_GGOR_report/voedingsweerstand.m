b = 200;
c = 100;
D = 10;
k = 10;
kD = k*D;
lam = sqrt(kD*c);
w = 0;

f = @(k,c) w/c * b/D + b./sqrt(k*D*c) .* coth( b./(k*D*c) ) - 2.0;

k = logspace(-5,5,11)

f(k,c)

%%

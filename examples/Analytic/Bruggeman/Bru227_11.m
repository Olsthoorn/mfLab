%% Bruggeman 227.11
% TO 121220

kD    = 500;
c     = 250;
S     = 0.001;
h     = 2;
R     = 300;
t     = 1;
beta2 = S/kD;
L2    = kD*c;
eta   = 1/(c*S);
r     = sinespace(R,5*R,50,pi/20,pi/2);

u     = logspace(-10,10,1000);
du    = diff(u); du = [0 du/2] + [du/2 0];
A     = (u.^2+R^2/L2).*bessely(0,u) - 2*S.*u.*bessely(1,u);
B     = (u.^2+R^2/L2).*besselj(0,u) - 2*S.*u.*besselj(1,u);

for ir = length(r):-1:1
    s     = (A.*besselj(0,r(ir)/R.*u)-B.*bessely(0,r(ir)/R.*u))./(A.^2+B.^2).*u;
    E     = exp(-eta*t-(u.^2*t)./(beta2*R^2));
    phi(ir) = 2*h/pi * sum(s.*E .*du);
end

figure; xlabel('r [m]'); ylabel('drawdown [m]'); title('Bruggeman 227.11');
plot(r,phi);
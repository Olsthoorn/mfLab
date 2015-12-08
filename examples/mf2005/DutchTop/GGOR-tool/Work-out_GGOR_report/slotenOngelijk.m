% GGOR report
% Nulpunt x-as in het midden van het perceel


L = 200; b = L/2;
kD =  50;
c  = 200;
lambda = sqrt(kD*c);

hL =  -1.0;
HR =   1.0;
phi = -1.0;

N  = 0.002;

x = -b:L/500:b;

xoL = x/lambda;
boL = b/lambda;

figure; set(gca,'nextplot','add');
xlabel('x [m]');
ylabel('elevation [m]');
title('Nulpunt x-as in het midden van perceel, sloten ongelijk, geen intredeweerstand')
for hR = hL:(HR-hL)/10:HR
    fprintf('hr = %g\n',hR);
    h = N*c + phi - (hL - hR)/2 * sinh(xoL)/sinh(boL) - (N*c+phi- (hL+hR)/2) * cosh(xoL)/cosh(boL);
    plot(x,h);
end

text(-0.8*b,0.5*(hL+hR),sprintf('kD = %.0f m^2/d\nc = %.0f d\nphi = %.2f m\nN = %.3fm/d',kD,c,phi,N));

kD = 4500;
S  = 0.10;
R0 = 4500;
t0 = R^2*S/(2.25 * kD);
Q1  = 2.5e6/365.24;
Q2  = 7.5e6/365.24;
r = logspace(2,log10(8500),50);
%figure('name','Theis','pos',screenPos(0.75));
hold on;

%t = t0:30:t0+180;
t = 0.61*365.25;
for it=1:numel(t)
    s = 1/(4*pi*kD) * expint(r.^2 * S / (4 * kD *t(it)));
    plot(r,Q1*s,'r--','lineWidth',3);
    plot(r,Q2*s,'b--','lineWidth',3);
end
%set(gca,'xScale','log','yDir','reverse');

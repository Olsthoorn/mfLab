R=1000; kD=1500;
x= [700,-500,900, 39,  -460,  -679,  -197,    -4,  -206];
y= [200, 750, 50, 583,    4,  -521,  -600,  -714,  -758];
Q=-[800,500,1200, 600, 1000,  300,   450,    450,   300];
kD=1500;
theta=[0:72]*2*pi/72; plot(R*cos(theta),R*sin(theta),'k','linewidth',2); hold on; %Plot the island
d1=sqrt(x.^2+y.^2); d2=R^2./d1; xm=x.*d2./d1; ym=y.*d2./d1; % d1, d2 and mirror wells
xg=-1500:10:1500; yg=-1500:10:1500; [X,Y]=meshgrid(xg,yg);  % grid points
%Superposition
s=zeros(size(X));
t=zeros(size(X));  % stroom functie
for i=1:length(x)
s=s+Q(i)/(2*pi*kD)*log(sqrt((X-x(i)).^2+(Y-y(i)).^2)./sqrt((X-xm(i)).^2+(Y-ym(i)).^2)*d1(i)/d2(i));
t=t+Q(i)*atan2(Y-y(i),X-x(i))-Q(i)*atan2(Y-ym(i),X-xm(i)); % doesn't work well
end
contour(X,Y,s,50); axis 'square'
contour(X,Y,t,50);

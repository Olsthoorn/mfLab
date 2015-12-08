function clrmap=jet3
%JET3 colormap with grey in the middle red on top and blue on bottom
%
% somewhat better suited for fresh-salt water in cross sections.
% used for several animations in SWIM conferences and on YouTube
%
% TO 111005

r=[0   1   1 0.5  ; 0/8 5/8 7/8  8/8];
g=[0   0 0.8 0.8 0 0; 0 1/8 3/8 5/8 7/8 8/8];
b=[0.5   1   1   0 ; 0 1/8 3/8 8/8];

% figure; hold on
% plot(r(2,:),r(1,:),'r');
% plot(g(2,:),g(1,:),'g');
% plot(b(2,:),b(1,:),'b');

L=64; x=0:L;x=0.5*(x(1:end-1)+x(2:end));

dR=diff(r,1,2);
red=zeros(1,L);
for i=1:size(dR,2)
    I=find(x>=L*r(2,i) & x<=L*r(2,i+1));
    red(I)=r(1,i)+dR(1,i)/dR(2,i)*(x(I)-x(I(1)))/L;
end        

dG=diff(g,1,2);
green=zeros(1,L);
for i=1:size(dG,2)
    I=find(x>=L*g(2,i) & x<=L*g(2,i+1));
    green(I)=g(1,i)+dG(1,i)/dG(2,i)*(x(I)-x(I(1)))/L;
end        

dB=diff(b,1,2);
blue=zeros(1,L);
for i=1:size(dB,2)
    I=find(x>=L*b(2,i) & x<=L*b(2,i+1));
    blue(I)=b(1,i)+dB(1,i)/dB(2,i)*(x(I)-x(I(1)))/L;
end        
% figure; hold on;
% plot(x,red  ,'r','Linewidth',2);
% plot(x,green,'g','Linewidth',2);
% plot(x,blue,'b' ,'Linewidth',2);

clrmap=[red' green' blue'];


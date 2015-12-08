% Check balance of model radmod1
% TO 051124

spd=86400; % seconds per day

close all

theta=1;

t=logspace(-2,4,31); % seconds

rw=0.02; zw=3.0;  % radius and depth of piezometer (flat open bottom)
rw=1;

r=logspace(-3,1,41);         rm=0.5*(r(1:end-1)+r(2:end));            % r from 1 mm to 10 m
z=flipud([-5:0.1:5]');       zm=0.5*(z(1:end-1)+z(2:end)); 

Nr=length(r); Nz=length(z); Nt=length(t);

KR=10/spd; SS=1e-5;
kr=KR*ones(Nz-1,Nr-1); kz=kr; % m/s
Ss=SS*ones(Nz-1,Nr-1);

FH=zeros(Nz,Nr)*NaN; FQ=zeros(Nz,Nr); Phi0=zeros(Nz,Nr);
%FH(:,end)=0;

% % wellpoint at -3 m alle cellen boven -3 en radius < 2 cm zijn dicht
Ir=find(rm<=rw);
Jw=find(zm>=zw);

%kr(Jw,Ir)=KR*1e6; kz(Jw,Ir)=KR*1e6; kr(Jw,Ir(end)+1)=0, kz(Jw,Ir(end)+1)=0;
%Ss(Jw(1),Ir)=100; 

J=find(z>=zw); I=find(r<=rw); Phi0(:,I)=1;

%FH(Jw(1),I)=-1;
%FQ(1,I(end-1))=10/spd;  % m3/s

[Phi,Q]=radmodt(r,z,t,kr,kz,Ss,FH,FQ,Phi0,theta);

P=permute(Phi,[3,2,1]);
for i=1:Nr
    semilogx(t,P(:,i,j)); hold on;
end

%figure
%PhiRange=[min(Phi(:)):(max(Phi(:))-min(Phi(:)))/30:max(Phi(:))];
%contour(r,z,Phi(:,:,end),PhiRange);
min(Phi(:)),max(Phi(:))
%set(gca,'xlim',[0,50],'ylim',[0,50]);

figure
for k=1:length(t)
plot(z,Phi(:,1,k)); hold on
end


pi*rw^2*max(Ss(:))*(Phi(1,1,1)-Phi(1,1,end))
pi*r(end)^2*abs(z(end)-z(1))*min(Ss(:))*(Phi(end,end,end)-Phi(end,end,1))

pi*rw^2*max(Ss(:))*Phi(1,1,1)  +pi*(r(end)^2-rw^2)*min(Ss(:))*Phi(end,end,1)
pi*rw^2*max(Ss(:))*Phi(1,1,end)+pi*(r(end)^2-rw^2)*min(Ss(:))*Phi(end,end,end)
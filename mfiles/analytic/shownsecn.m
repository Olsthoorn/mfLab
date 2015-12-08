function shownsecn(phi,q,s,Z,h,X,x)
%SHOWNSECN shows results of nsecn
%
% Example:
%    shownsecn(phi,q,s,Z,h,X,x) 
%
% Results are in totals per aquifer ans as a spatial contour plot of head and stream lines
% Based on the stream function.
%
% See also: nsecn hantushn
%
% TO 000112  001202


if size(Z,1)<2*size(phi,1)-1, 
    error('All layer interfaces must be specified length(z(:,1)) is %d but must be %d',size(Z,1),size(phi,1)*2+1);
end

nLay=size(phi,1);
linestyle={'-','--',':','-.'};
linecolor={'b','r','g','k','m','c','y'};
H=zeros(size(X));
for i=1:length(x)-1,
   I=find(X>=x(i) & X<=x(i+1));
   if ~isempty(I),
      H(I)=h(i);
   end
end

dd=[0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000 2000 5000 10000];

figure;
a1=subplot(3,1,1); hold on; grid on; ylabel('[m]');    title('heads in all aquifers');
a2=subplot(3,1,2); hold on; grid on; ylabel('[m2/d]'); title('flux, positive = --> (black is total q)');
a3=subplot(3,1,3); hold on; grid on; ylabel('[m/d]');  title('seepage from top'); xlabel('x [m]');

figure;
b1=axes; hold on; xlabel('x [m]'); ylabel('z [m]'); title('head- and streamlines');


axes(a1); hH=plot(X,H,'k');        legh=[]; hh=plot(X,phi);
axes(a2); hQ=plot(X,sum(q,1),'k'); legq=[]; hq=plot(X,q);
axes(a3);                          legs=[]; hs=plot(X,s);
for i=1:nLay
   j=rem(i,length(linestyle)); if j==0, j=length(linestyle); end
   k=rem(i,length(linecolor)); if k==0, k=length(linecolor); end
   set(hh(i),'linestyle',linestyle{j},'color',linecolor{k});
   set(hq(i),'linestyle',linestyle{j},'color',linecolor{k});
   set(hs(i),'linestyle',linestyle{j},'color',linecolor{k});
   legh{i}=['h',sprintf('%d',i)];
   legq{i}=['q',sprintf('%d',i)];
   legs{i}=['s',sprintf('%d',i)];
end
axes(a1); legend([{'hT'},legh]);
axes(a2); legend([{'qT'},legq]);
axes(a3); legend(legs);

axes(b1);
mPsi=min(0,min(sum(q,1))); MPsi=max(0,max(sum(q,1)));
mPhi=min([phi(:);h(:)]);   MPhi=max([phi(:);h(:)]);
mPsi=floor(mPsi); MPsi=ceil(MPsi); dPsi=(MPsi-mPsi)/20; dPsi=dd(max(find(dPsi>=dd))); psirange=mPsi:dPsi:MPsi;
mPhi=floor(mPhi); MPhi=ceil(MPhi); dPhi=(MPhi-mPhi)/50; dPhi=dd(max(find(dPhi>=dd))); phirange=mPhi:dPhi:MPhi;


if size(Z,2)==1,    Z=Z*ones(size(x));   end

FAR=1e6;
Z=interp1([X(1)-FAR,x(2:end-1),X(end)+FAR]',Z',X')';

Z(1,:)=H;


cont(phi,q,ones(size(Z(:,1)))*X,Z,H,psirange,phirange);
axes(b1);
set(b1,'xLim',[min(X(:)),max(X(:))],'ylim',[min(Z(:)),max(Z(:))]);
title(sprintf('head- and streamlines, dPhi=%.2f m, dPsi=%.2f m2/d',dPhi,dPsi));

function cont(phi,q,x,z,H,psirange,phirange)
N=length(phi(:,1));
Phi=zeros(2*N+1,length(x));
Psi=zeros(2*N+1,length(x));

Q=flipud(cumsum(flipud(q),1));
Psi(1:2:2*N-1,:)= Q;
Psi(2:2:2*N,  :)= Q;
Psi(2*N+1,    :)=zeros(size(Q(end,:)));
Phi(1,        :)= H;
Phi(2:2:2*N,  :)= phi;
Phi(3:2:2*N+1,:)= phi;

for i=1:length(z(:,1))-1
   if rem(i,2)
      if i==1
          patch([x(i,1:end),fliplr(x(i+1,1:end))],...
              [H(i,1:end),fliplr(z(i+1,1:end))],[0.85,0.85,0.85]);
      else
          patch([x(i,1:end),fliplr(x(i+1,1:end))],...
              [z(i,1:end),fliplr(z(i+1,1:end))],[0.85,0.85,0.85]);
      end
   else
      patch([x(i,1:end),fliplr(x(i+1,1:end))],...
         [z(i,1:end),fliplr(z(i+1,1:end))],[1,1,0.5]);
   end
end

z=min(Phi,z);

contour(x(1:end,:),z(1:end,:),Phi(1:end,:),phirange,'r');
contour(x(1:end,:),z(1:end,:),Psi(1:end,:),psirange,'b');
line(x,Phi(3,:),'color','r','linewidth',1);
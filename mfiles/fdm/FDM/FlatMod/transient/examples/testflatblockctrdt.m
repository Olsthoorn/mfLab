%test flatblockctrd1
% TO 010828

clear; close all
k=2;
switch k
case 1  % steady state analytic solution
   x=[-0.1,0:100:2000,2000.1]; y=[0,100];
   n=0.01;
case 2  % flow to a well, theis test
   x=[0,logspace(0,4,40)];
   y=[0,logspace(0,4,40)];
	n=0.0;	%10 mm/d;
case 3
end
   
% General: compute centers of cells   
xm=0.5*(x(1:end-1)+x(2:end)); Nx=length(x); dx=diff(x);
ym=0.5*(y(1:end-1)+y(2:end)); Ny=length(y); dy=abs(diff(y(:)));

% time
t=[0,logspace(0,4,40)]; tm=0.5*(t(1:end-1)+t(2:end)); Nt=length(t); dt=diff(t);

% conductivities and storage (transmissivities and storage)
kx=100*ones(Ny-1,Nx-1);
ky=kx;
S =0.1*ones(Ny-1,Nx-1);

S(10:15,10:15)=NaN;						% this is to make some cells inactive

FQ=dy*dx*n; FQ(1,1)=-600;				% fixed Q, (one quarter model, so copare wirth 2400 m3/d)

FH=NaN*FQ; %FH(1)=0;						% fixed head boundaries

IH=FH; IH(find(isnan(FH)))=0;			% initial heads

GHB=[];										% not used now, use GHB=dy*dx./c if c is the top system resistance

theta=0.75;									% implicitness (0-1)
solver=1;									% 1 is direct, 2 is pcg with preconditioning (may be faster for large models)

[Phi,Q]=flatblockctrdt(x,y,t,kx,ky,S,FH,FQ,GHB,IH,theta,solver);	% run the model

% RESULTS
switch k
case 1		% analytic 1 d strip, converges to final solution, ok
	figure,
	for i=1:length(t)
   	h=Phi(1,:,i); plot(xm,h); hold on
   end
   L=4000; yy=x.*(L-x)*n./(2*kx(1));   
   plot(x,yy,'r');
   xlabel('x [m]'); ylabel('head or drawdown [m]'); title('strip of land');
case 2		% theis test, ok (note the line above, allowing to make cells inactive)
   
   % comparing with Theis solution
   mu=S(1,1); kr=kx(1); D=1; r=x(2:end); tt=t(2:end);							% take data from numeric data
   u=mu/(4*kr*D)*(ones(size(tt(:)))*r.^2)./(tt(:)*ones(size(r)));			% argument of Theis function
	phiTh=4*FQ(1,1)./(4*pi*kr*D)*expint(u);										% theis solution
   plot(r,phiTh','r');																	% plot Theis lines for each time
   
   % show numeric soltion
   for i=1:length(t)
	   h=Phi(1,:,i); hold on															% get head along xaxis for this time
      plot(xm,h);
   end
   set(gca,'xscale','log');															% semilog for better visibility
   xlabel('x [m]'); ylabel('head [m]'); title('Theis, drawdown versus time (red analytic, blue model)');
   % contour the heads of the end time
   figure
   contour(xm,ym,Phi(:,:,end),[-40:0.5:0]); xlabel('x [m]'), ylabel('y [m]'); title('Theis')
end

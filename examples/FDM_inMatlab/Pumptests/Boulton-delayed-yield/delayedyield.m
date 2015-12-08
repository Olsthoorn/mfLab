% dimensionless parameters
% The idea was to compute type curves for delayed yield by means of a
% numerical axisymmetric model. However, as it turns out, the values W(u) thus
% computed that are lower than about 0.01 are inaccurate, probably because
% of numerical diffictulties, even with Matlab's backslash solution.
% This makes this
% method highly unrealiable, even though it may not matter too much in
% practice.

% A better way might be using a multi-layer analytical Hantush function,
% under the condition that the necessary back transformation from the
% Laplace space is accurate. In that case, we have are free to establish
% a full 3d multilayer model as was developed by Kees Maas and Kick Hemker
%
% The model an be used to obtain also delayed yield by having an upper
% layer to provide the delayed water from its storage and a second layer to
% provide the specific yield as such.

close all
clear variables

% Theis curve W(u)
% s/(Q/(4*pi*T)) = W(u), u=r^2*S/(4Tt);

% variables, which can be arbitrry but we use them such that
% Q/(4piT)=1 and S/(4T)=1;
T=1;        % arbirtrary
Q=4*pi*T;   % so that s/(Q/(4piT))==s=W
S=4*T;      % so that S/4T ==1
r0=1;       % so that r^2S/(4Tt)=u=1/t

r=logspace(-2,5,71);
t=logspace(-3,3,61)';

y=[-1 0];
[r,y,rm,ym,dr,dy,Nr,Nm]=modelsize(r,y);

u=t*(1./rm).^2;

%IH=zeros(size(rm));

u=S./(4*T*t)*(rm.^2);
W=expint(u);
s=Q/(4*pi*T)*W;

IH=s(1,:); % proper analytical start values for the numerial model
FH=NaN(size(rm));
FQ=zeros(size(rm)); FQ(1)=Q;

[Phi,Qt,Qx,Qy,Qs]=fdm2t(r,y,t,T,T,S,IH,FH,FQ,'axisymmetric'); Phi=permute(Phi,[3,2,1]);

Phi(Phi<1e-10)=NaN;
s  (s  <1e-10)=NaN;

figure; hold on; set(gca,'xscale','log','yscale','log');
plot(t,s,'b');
plot(t,Phi,'r');

%%
figure; hold on; set(gca,'xscale','log','yscale','log');
plot(t,s,'b');
plot(t,Phi,'r');
plot(1./u,W,'go');

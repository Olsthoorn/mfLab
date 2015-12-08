% Some multilayer solutions of Bruggeman (1999) solutions in series 710
% TO 101204. Adapted to Octave (3.6.2) PRN 111105

clear variables;
close all;

global transmissivity

isOctave = exist ('OCTAVE_VERSION', 'builtin') ~= 0;
isOct = false;
if (isOctave)
  OV = strsplit (OCTAVE_VERSION, '.');
  isOct = (str2num (OV{2}) >= 6);
end

coshm=inline('(expm(x)+expm(-x))/2');
sinhm=inline('(expm(x)-expm(-x))/2');

%% given values

b=50;
NLay=5;  % choose n layers

z0 = 0;  % top of aquifer aquitard system
d  = [   5   10    5  10   2   10    5    2     5   10]; % aquitard thickness
D  = [  15   20   15  20  20   15   10   10    15   20]; % aquifer thickness
q  = [   5   -2    1  -1   3   10   -5    7     1   -6]*0.001; % discharge
Q  = [   1    2    0  -3   1   -2    1    2   -10    3]; % discharge flat
h  = [   1   -1    0  -3   1   -2    1    2   -10    3]; % extraction radial
k  = [   2    2    5  10  15    1    5   10   0.5   10]; % conductivity
c  = [100    10  200  40 100    5   10    3    10  100]; % resistance confining beds
w  = [   1    5    1   2   2    3    1    1     2    1]; % entry resistance
S  = [   1    1    1   1   1    1    1    1     1    1]*1e-3; % elastic storage

q=abs(q);
Q=abs(Q);
h=abs(h);

kD = k.*D;         % transmissivity

DZ=reshape([d;D],[2*length(d),1]); Z=z0-[0; cumsum(DZ)];
Z=Z(1:2*NLay+1);
Iz=sort([1:NLay 1:NLay]);
%% Derive matrices

T  = diag(kD(1:NLay));   T_m1=T^(-1);  % trasmissivity matrix
S_m1= diag(1./S);     % storage ceofficient matrix
H   = diag(1./(k(1:NLay).*w(1:NLay))); H_m1= H^(-1); % entry resistance matrix
I  = eye(NLay);

A=-diag( 1./(kD(2:NLay  ).*c(2:NLay)),-1)+...
  +diag( 1./(kD(1:NLay  ).*c(1:NLay))+1./(kD(1:NLay).*[c(2:NLay) Inf]), 0)+...
  -diag( 1./(kD(1:NLay-1).*c(2:NLay)),1);

sqrtA  = sqrtm(A);

ATm1q=A\(T_m1*q(1:NLay)');              % simple vector in many solutions
sqrtATm1Q=(sqrtA\(T_m1*Q(1:NLay)'));    % simple vector in many solutions
sqrtATm1q=(sqrtA\(T_m1*q(1:NLay)'));    % same used in axial symmetric solutions (720 series)

snhm=sinhm(b*sqrtA); snhm_m1=snhm^(-1);
cshm=coshm(b*sqrtA); coshm_m1=cshm^(-1);

cS=diag(c(1:NLay).*S(1:NLay)); cS_m1=cS^(-1);

x=0:b;               Nx=length(x);
t=logspace(-2,1,41); Nt=length(t);
Phit=zeros(NLay,Nt);
Qt  =zeros(NLay,Nt);
Phix=zeros(NLay,Nx);
Qx  =zeros(NLay,Nx);


transmissivity=kD(1:NLay)';

%% Solution 710.01 =====================================================

% voldoet aan de differentiaalvergelijking, maar de dimensie van A is hier
% 1/d want de vector h heeft elementen h/(c S).

s={ 'Bruggeman (1999) Solution 710.01. Transient.';
    'n-layer system under open water with sudden rise h of its level';
    'phi=phi(t). Only vertical flow in semi-pervious layers.'
    };
if (~isOct); s = s{1}; end

AA=-diag( 1./(S(2:NLay  ).*c(2:NLay)),-1)+...
   +diag( 1./(S(1:NLay  ).*c(1:NLay))+1./(S(1:NLay).*c(2:NLay+1)), 0)+...
   -diag( 1./(S(1:NLay-1).*c(2:NLay)),1);

AA_m1=AA^(-1);

htop=h; htop(2:end)=0;

for it=1:length(t)
    Phit(:,it)=(I-expm(-AA*t(it)))*AA_m1*cS_m1*htop(1:NLay)';
    % Qt remains every where zeros, there is no horizontal discharge;
end

figure;
subplot(2,1,1); hold on; grid on; plot(t,Phit);
title(s); xlabel('time [d]'); ylabel('drawdown [m]');

%% Solution 710.02 =====================================================

s={ 'Bruggeman (1999) Solution 710.01. Steady state.';
    'n-layer system under open water with sudden rise h of its level';
    'phi=phi(t). Only vertical flow in semi-pervious layers.'
    };
if (~isOct); s = s{1}; end

htop=h; htop(2:end)=0;

for it=1:length(t)
    Phit=AA_m1*cS_m1*htop(1:NLay)';
end

plot(t([1 end]),[Phit Phit]); grid on; title(s);

%% Solution 710.12 =====================================================

s={ 'Bruggeman (1999) Solution 710.02. Steady state.';
    'All aquifers with open boundary. Sudden drawdown of the surface water level,';
    'which is kept constant thereafter. phi=phi(x)=drawdown'
    };

x=0:b; Nx=length(x);

Phix=zeros(NLay,Nx);
Qx  =zeros(NLay,Nx);

for ix=1:length(x)
    Phix(:,ix)=  expm(-x(ix)*sqrtA)*      h(1:NLay)';
    Qx(:,ix)  =T*expm(-x(ix)*sqrtA)*sqrtA*h(1:NLay)';
end

plotbrug(x,Z,Phix,Qx,s);

%% Solution 710.12 =====================================================

s={ 'Bruggeman (1999) Solution 710.12. Steady state.';
    'All aquifers with open boundary. Sudden drawdown of the surface water level,';
    'which is kept constant thereafter. phi=phi(x)=drawdown'
    };

x=0:b; Nx=length(x);

Phix=zeros(NLay,length(x));
Qx  =zeros(NLay,length(x));
for ix=1:length(x)
    Phix(:,ix)=  expm(-x(ix)*sqrtA)*h(1:NLay)';
    Qx(  :,ix)=T*expm(-x(ix)*sqrtA)*sqrtA*h(1:NLay)';
end

plotbrug(x,Z,Phix,Qx,s);

%% Solution 710.13 =====================================================

s={ 'Bruggeman (1999) Solution 710.13. Steady state.';
    'All aquifers with open boundary. The plane at x=0 ketp at zero head.';
    'Constant vertical infiltration q into the aquifers.'
    }; % phi=phi(x)=drawdown';

x=logspace(0,3,31); Nx=length(x);

Phix=zeros(NLay,Nx);
Qx  =zeros(NLay,Nx);
for ix=1:length(x)
    Phix(:,ix)=(I-expm(-x(ix)*sqrtA))*ATm1q;
    Qx(  :,ix)=-T*expm(-x(ix)*sqrtA)*sqrtA*ATm1q;
end

plotbrug(x,Z,Phix,Qx,s);

%% Solution 710.15 =====================================================

s= { 'Bruggeman (1999) Solution 710.14. Steady state.';
     'Fully penetrating well screens in all aquifers at x=0.';
     'Sudden discharges which are kept constant theraafter.'
    }; % phi=phi(x)=drawdown';

x=logspace(0,3,31); Nx=length(x);

Phix=zeros(NLay,Nx);
Qx  =zeros(NLay,Nx);
for ix=1:Nx
    Phix(:,ix)=  expm(-x(ix)*sqrtA)*      sqrtATm1Q/2;
    Qx(:,ix) = T*expm(-x(ix)*sqrtA)*sqrtA*sqrtATm1Q/2;
end

plotbrug(x,Z,Phix,Qx,s);

%% Solution 710.16 =====================================================

s={ 'Bruggeman (1999) Solution 710.16. Steady state.';
    'All aquifers with open boundary with entrance resistance.';
    'Sudden drawdown of the surface water level which i kept constant thereafter.'
    }; % phi=phi(x)=drawdown';

x=logspace(0,3,31); Nx=length(x);

Phix=zeros(NLay,length(x));
Qx  =zeros(NLay,length(x));

sqrtApHm1Hh=(sqrtA+H)\(H*h(1:NLay)');

for ix=1:length(x)
    Phix(:,ix)=  expm(-x(ix)*sqrtA)*      sqrtApHm1Hh;
    Qx(:,ix)  =T*expm(-x(ix)*sqrtA)*sqrtA*sqrtApHm1Hh;
end

plotbrug(x,Z,Phix,Qx,s);

%% Solution 710.17 =====================================================

s={ 'Bruggeman (1999) Solution 710.17. Steady state.';
    'All aquifers with open boundary with entrance resistance.';
    'Constant infiltration. Zero head at x=0.'
    }; % phi=phi(x)=drawdown';

x=logspace(0,3,31); Nx=length(x);

Phix=zeros(NLay,length(x));
Qx  =zeros(NLay,length(x));

B=(sqrtA+H)\(H*ATm1q);
for ix=1:length(x)
    Phix(:,ix)= ATm1q-expm(-x(ix)*sqrtA)*B;
    Qx(  :,ix)=    -T*expm(-x(ix)*sqrtA)*sqrtA*B;
end

plotbrug(x,Z,Phix,Qx,s);

%% Solution 710.18 =====================================================

s={ 'Bruggeman (1999) Solution 710.18. Steady state.';
    'Infiltration into aquifers with q(x)=q for |x|<b en 0 for |x|>b.';
    'Flux=0 for x=0.'
    }; % phi=phi(x)=drawdown';

x=unique([0 b logspace(0,log10(3*b),41)]); x=unique([-x x]);

bexpm  = expm(-b*sqrtA);

x=0:100; Nx=length(x);

Phix=zeros(NLay,length(x));
Qx  =zeros(NLay,length(x));

for ix=1:length(x)
    if abs(x(ix))<=b
        Phix(:,ix)=(I-coshm(x(ix)*sqrtA)*     bexpm)*ATm1q;
        Qx(  :,ix)=T*sinhm(x(ix)*sqrtA)*sqrtA*bexpm *ATm1q;
    else
        Phix(:,ix)=  expm(-abs(x(ix))*sqrtA)*      snhm*ATm1q;
        Qx(  :,ix)=T*expm(-abs(x(ix))*sqrtA)*sqrtA*snhm*ATm1q;
    end
end

plotbrug(x,Z,Phix,Qx,s);

%% Solution 710.21 =====================================================

s={ 'Bruggeman (1999) Solution 710.25. Steady state.';
    'Given drawdown h for x=b. Flux=0 for x=0.'
     }; % phi=phi(x)=drawdown';

x=0:100; Nx=length(x);

Phix=zeros(NLay,length(x));
Qx  =zeros(NLay,length(x));

for ix=1:length(x)
    Phix(:,ix)=   coshm(x(ix)*sqrtA)*      coshm_m1*h(1:NLay)';
    Qx(  :,ix)=-T*sinhm(x(ix)*sqrtA)*sqrtA*coshm_m1*h(1:NLay)';
end

plotbrug(x,Z,Phix,Qx,s);

%% Solution 710.22 =====================================================

s={ 'Bruggeman (1999) Solution 710.22. Steady state.';
    'Given drawdown h for x=b. Zero drawdown at x=0.'
     }; % phi=phi(x)=drawdown';

x=0:b; Nx=length(x);

Phix=zeros(NLay,length(x));
Qx  =zeros(NLay,length(x));

for ix=1:length(x)
    Phix(:,ix)=   sinhm(x(ix)*sqrtA)*      snhm_m1*h(1:NLay)';
    Qx(  :,ix)=-T*coshm(x(ix)*sqrtA)*sqrtA*snhm_m1*h(1:NLay)';
end

plotbrug(x,Z,Phix,Qx,s);

%% Solution 710.23 =====================================================

s={ 'Bruggeman (1999) Solution 710.23. Steady state.';
    'Zero head at x=b. Zero flux at x=0. Constant infitlration.'
     }; % phi=phi(x)=drawdown';

x=0:b; Nx=length(x);

Phix=zeros(NLay,length(x));
Qx  =zeros(NLay,length(x));

for ix=1:length(x)
    Phix(:,ix)=(I-coshm(x(ix)*sqrtA)*      coshm_m1)*ATm1q;
    Qx(  :,ix)= T*sinhm(x(ix)*sqrtA)*sqrtA*coshm_m1 *ATm1q;
end

plotbrug(x,Z,Phix,Qx,s);

%% Solution 710.24 =====================================================

s={ 'Bruggeman (1999) Solution 710.24. Steady state.';
    'Fully penetrating well screens in all aquifers at x=0, -b, 3b, -3b etc.';
    'Constant but different discharges q. Flux=0 for x=0.'
     }; % phi=phi(x)=drawdown';

x=0:b; Nx=length(x);

Phix=zeros(NLay,length(x));
Qx  =zeros(NLay,length(x));
for ix=1:length(x)
    Phix(:,ix)=   coshm(x(ix)*sqrtA)*      snhm_m1*sqrtATm1Q/2;
    Qx(  :,ix)=-T*sinhm(x(ix)*sqrtA)*sqrtA*snhm_m1*sqrtATm1Q/2;
end

plotbrug(x,Z,Phix,Qx,s);

%% Solution 710.25 =====================================================

s={ 'Bruggeman (1999) Solution 710.25. Steady state.';
    'Fully penetrating well screen in all aquifers in the middle of a strip with width 2b.';
    'Constant but different dischargfes q. Drawdown=0 for x=0 and x=2b.'
    }; % phi=phi(x)=drawdown';

x=0:b; Nx=length(x);

Phix=zeros(NLay,length(x));
Qx  =zeros(NLay,length(x));

for ix=1:length(x)
    Phix(:,ix)=   sinhm(x(ix)*sqrtA)*      coshm_m1*sqrtATm1Q/2;
    Qx(  :,ix)=-T*coshm(x(ix)*sqrtA)*sqrtA*coshm_m1*sqrtATm1Q/2;
end

plotbrug(x,Z,Phix,Qx,s);

%% Solution 710.26 =====================================================

s={ 'Bruggeman (1999) Solution 710.25. Steady state. ';
    'All aquifers with entrance resistances w_i to open water with a constant drawdown.';
    'Zero flux at x=0. phi=phi(x)=drawdown.'
    };

F_m1=(H_m1*sqrtA*snhm+cshm)^(-1);

x=0:b; Nx=length(x);

Phix=zeros(NLay,length(x));
Qx  =zeros(NLay,length(x));
for ix=1:length(x)
    Phix(:,ix)=   coshm(x(ix)*sqrtA)*      F_m1*h(1:NLay)';
    Qx(  :,ix)=-T*sinhm(x(ix)*sqrtA)*sqrtA*F_m1*h(1:NLay)';
end

plotbrug(x,Z,Phix,Qx,s);

%% Solution 710.27 =====================================================

s={ 'Bruggeman (1999) Solution 710.27. Steady state.';
    'All aquifers with entrance resistance w_i to oen water at x=b.';
    'Constant drawdown of the open water level. Zero drawdown at x=0.'
    }; % phi=phi(x)=drawdown';

F_m1=(H_m1*sqrtA*cshm+snhm)^(-1);

x=0:b; Nx=length(x);

Phix=zeros(NLay,length(x));
Qx  =zeros(NLay,length(x));
for ix=1:length(x)
    Phix(:,ix)=   sinhm(x(ix)*sqrtA)*      F_m1*h(1:NLay)';
    Qx(  :,ix)=-T*coshm(x(ix)*sqrtA)*sqrtA*F_m1*h(1:NLay)';
end

plotbrug(x,Z,Phix,Qx,s);

%% Solution 710.28 =====================================================

s={
    'Bruggeman (1999) solution 710.28. Steady state.';
    'Alle aquifers with entrance resistance w to open water at x=b.';
    'Vertical infiiltration q. Flux=0 for x=0.'
    };  % phi=phi(x)=drawdown';

F_m1=(H_m1*sqrtA*snhm+cshm)^(-1);

x=0:b; Nx=length(x);

Phix=zeros(NLay,length(x));
Qx  =zeros(NLay,length(x));
for ix=1:length(x)
    Phix(:,ix)=(I-coshm(x(ix)*sqrtA)*      F_m1)*ATm1q;
    Qx(  :,ix)= T*sinhm(x(ix)*sqrtA)*sqrtA*F_m1 *ATm1q;
end

plotbrug(x,Z,Phix,Qx,s);

%% Solution 720.01 =====================================================

R=b;

[V,D]=eig(R*sqrtA); V_m1=V^(-1);
K0R=V*diag(besselk(0,diag(D)))*V_m1; K0R_m1=K0R^(-1);
K1R=V*diag(besselk(1,diag(D)))*V_m1; K1R_m1=K1R^(-1);
I0R=V*diag(besseli(0,diag(D)))*V_m1; I0R_m1=I0R^(-1);
I1R=V*diag(besseli(1,diag(D)))*V_m1; I1R_m1=I1R^(-1);

s={
    'Bruggeman (1999) solution 720.01. Steady state.';
    'Vertical infiiltration q(r). q(r)=q for 0<r<=R and 0 for r>R.'
    };  % phi=phi(x)=drawdown';

r=unique([b logspace(0,log10(5*b),31)]); Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);

for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    if r(ir)<=R
        Phir(:,ir)=(I-          R*sqrtA*K1R*v*diag(besseli(0,d))*v_m1)*      ATm1q;
        Qr(  :,ir)=2*pi*r(ir)*T*R*sqrtA*K1R*v*diag(besseli(1,d))*v_m1*sqrtA *ATm1q;
    else
        Phir(:,ir)=             R*I1R*v*diag(besselk(0,d))*v_m1*      sqrtATm1q;
        Qr(  :,ir)=2*pi*r(ir)*T*R*I1R*v*diag(besselk(1,d))*v_m1*sqrtA*sqrtATm1q;
    end
end

plotbrug(r,Z,Phir,Qr,s);


%% Solution 720.03 =====================================================

s={
    'Bruggeman (1999) solution 720.03. Steady state.';
    'Fully penetrating wells in al aquifers at r=0.';
    'Constant but different discharges Q'
    };  % phi=phi(x)=drawdown';

r=unique([b logspace(0,log10(5*b),31)]); Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);

for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    Phir(:,ir)=             v*diag(besselk(0,d))*v_m1*      T_m1*q(1:NLay)';
    Qr(  :,ir)=2*pi*r(ir)*T*v*diag(besselk(1,d))*v_m1*sqrtA*T_m1*q(1:NLay)';
end

plotbrug(r,Z,Phir,Qr,s);

%% Solution 720.11 =====================================================

s={
    'Bruggeman (1999) solution 720.11. Steady state.';
    'All aquifers with open boundary. Constant drawdown of the open water level.'
    };  % phi=phi(x)=drawdown';

r=unique([b logspace(0,log10(5*b),31)]); Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);
for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    Phir(:,ir)=             v*diag(besselk(0,d))*v_m1*      K0R_m1*h(1:NLay)';
    Qr(  :,ir)=2*pi*r(ir)*T*v*diag(besselk(0,d))*v_m1*sqrtA*K0R_m1*h(1:NLay)';
end

plotbrug(r,Z,Phir,Qr,s);

%% Solution 720.12 =====================================================

s={
    'Bruggeman (1999) solution 720.12. Steady state.';
    'Constant lowering of the polder level around a circular basin.';
    'Zero drawdown at r=R.'
    };  % phi=phi(x)=drawdown';

r=unique([b logspace(0,log10(5*b),31)]); Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);

for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    Phir(:,ir)=           (I-K0R_m1*v*diag(besselk(0,d))*v_m1)*      h(1:NLay)';
    Qr(  :,ir)=-2*pi*r(ir)*T*K0R_m1*v*diag(besselk(1,d))*v_m1 *sqrtA*h(1:NLay)';
end

plotbrug(r,Z,Phir,Qr,s);

%========================================================================
% Octave kan niet meer dan zo'n 20 complexe plots tegelijk aan. Wis de eerste serie
if (isOctave)
  printf ('Eerste 20 plots gereed.\n');
  printf ('Na een druk op <ENTER> toets worden ze gewist en komt\n');
  printf ('de volgende serie...\n');
  pause
  close all
end

%% Solution 720.13 =====================================================

s={
    'Bruggeman (1999) solution 720.12. Steady state.';
    'All aquifers with open boundary and zero head at r=R. Constant infiltration q.'
    };  % phi=phi(x)=drawdown';

r=unique([b logspace(0,log10(5*b),31)]); Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);

for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    Phir(:,ir)=           (I-K0R_m1*v*diag(besselk(0,d))*v_m1)*      ATm1q;
    Qr(  :,ir)=-2*pi*r(ir)*T*K0R_m1*v*diag(besselk(1,d))*v_m1 *sqrtA*ATm1q;
end

plotbrug(r,Z,Phir,Qr,s);

%% Solution 720.13 =====================================================

s={
    'Bruggeman (1999) solution 720.13. Steady state.';
    'Fully penetrating circular well screens with unilateral discharges q in all';
    'aquifers at r=R.'
    };  % phi=phi(x)=drawdown';

r=unique([b logspace(0,log10(5*b),31)]); Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);

for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    Phir(:,ir)=             v*diag(besselk(0,d))*v_m1*      K1R_m1*sqrtATm1Q;
    Qr(  :,ir)=2*pi*r(ir)*T*v*diag(besselk(0,d))*v_m1*sqrtA*K1R_m1*sqrtATm1Q;
end

plotbrug(r,Z,Phir,Qr,s);

%% Solution 720.15 =====================================================

s={
    'Bruggeman (1999) solution 720.15. Steady state.';
    'All aquifers with entrance resistance w to open water with a constant drawdown h at r=R.';
    'aquifers at r=R.'
    };  % phi=phi(x)=drawdown';

r=unique([b logspace(0,log10(5*b),31)]); Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);

F_m1=(H_m1*sqrtA*K1R+K0R)^(-1);
for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    Phir(:,ir)=             v*diag(besselk(0,d))*v_m1*      F_m1*h(1:NLay)';
    Qr(  :,ir)=2*pi*r(ir)*T*v*diag(besselk(1,d))*v_m1*sqrtA*F_m1*h(1:NLay)';
end

plotbrug(r,Z,Phir,Qr,s);

%% Solution 720.16 =====================================================

s={
    'Bruggeman (1999) solution 720.16. Steady state.';
    'All aquifers with entrance resistance w to open water with zero drawdown at r=R';
    'Constant lowering of the polder level around a circular basin.'
    };  % phi=phi(x)=drawdown';

r=unique([b logspace(0,log10(5*b),31)]); Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);

F_m1=(H_m1*sqrtA*K1R+K0R)^(-1);
for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    Phir(:,ir)=           (I-v*diag(besselk(0,d))*v_m1*      F_m1)*h(1:NLay)';
    Qr(  :,ir)=-2*pi*r(ir)*T*v*diag(besselk(1,d))*v_m1*sqrtA*F_m1 *h(1:NLay)';
end

plotbrug(r,Z,Phir,Qr,s);

%% Solution 720.17 =====================================================

s={
    'Bruggeman (1999) solution 720.17. Steady state.';
    'All aquifers with entrance resistance w to open water with zero drawdown at r=R';
    'Constant infiltration q.'
    };  % phi=phi(x)=drawdown';

r=unique([b logspace(0,log10(5*b),31)]); Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);

F_m1=(H_m1*sqrtA*K1R+K0R)^(-1);
for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    Phir(:,ir)=           (I-v*diag(besselk(0,d))*v_m1*      F_m1)*ATm1q;
    Qr(  :,ir)=-2*pi*r(ir)*T*v*diag(besselk(1,d))*v_m1*sqrtA*F_m1 *ATm1q;
end

plotbrug(r,Z,Phir,Qr,s);

%% Solution 720.21 =====================================================

s={
    'Bruggeman (1999) solution 720.21. Steady state.';
    'Given drawdown h for r=R. All aquifers with open boundary to the surface water.';
    'Flux=0 at r=0.'
    };  % phi=phi(x)=drawdown';

r=1:R; Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);
for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    Phir(:,ir)=              v*diag(besseli(0,d))*v_m1*      I0R_m1*h(1:NLay)';
    Qr(  :,ir)=-2*pi*r(ir)*T*v*diag(besseli(1,d))*v_m1*sqrtA*I0R_m1*h(1:NLay)';
end

plotbrug(r,Z,Phir,Qr,s);

%% Solution 720.22 =====================================================

s={
    'Bruggeman (1999) solution 720.22. Steady state.';
    'All aquifers with open boundary to the surface water with zero head at r=R.';
    'Flux=0 at r=0. Constant infiltration `.'
    };  % phi=phi(x)=drawdown';

r=1:R; Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);

for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    Phir(:,ir)=           (I-v*diag(besseli(0,d))*v_m1*      I0R_m1)*ATm1q;
    Qr(  :,ir)=2*pi*r(ir)*T*v*diag(besseli(1,d))*v_m1*sqrtA*I0R_m1 *ATm1q;
end

plotbrug(r,Z,Phir,Qr,s);

%% Solution 720.23 =====================================================

s={
    'Bruggeman (1999) solution 720.23. Steady state.';
    'All aquifers with entrance resistance w to open water with constant drawdown.';
    'Zero flux at r=0.'
    };  % phi=phi(x)=drawdown';

r=1:R; Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);

F_m1=(H_m1*sqrtA*I1R+I0R)^(-1);
for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    Phir(:,ir)=              v*diag(besseli(0,d))*v_m1*      F_m1*h(1:NLay)';
    Qr(  :,ir)=-2*pi*r(ir)*T*v*diag(besseli(1,d))*v_m1*sqrtA*F_m1*h(1:NLay)';
end

plotbrug(r,Z,Phir,Qr,s);

%% Solution 720.24 =====================================================

s={
    'Bruggeman (1999) solution 720.24. Steady state.';
    'All aquifers with entrance resistance w to open water with constant drawdown.';
    'Zero head at r=R. Zero flux at r=0. Constant infiltration q.'
    };  % phi=phi(x)=drawdown';

r=1:R; Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);

F_m1=(H_m1*sqrtA*I1R+I0R)^(-1);
for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    Phir(:,ir)=          (I-v*diag(besseli(0,d))*v_m1*      F_m1)*ATm1q;
    Qr(  :,ir)=2*pi*r(ir)*T*v*diag(besseli(1,d))*v_m1*sqrtA*F_m1 *ATm1q;
end

plotbrug(r,Z,Phir,Qr,s);

%% Solution 720.31 =====================================================

s={
    'Bruggeman (1999) solution 720.31. Steady state.';
    'Fully penetrating wells in all aquifers at r=0.';
    'Constant but different discharges Q. Open water boundary at r=R with zero drawdown.'
    };  % phi=phi(x)=drawdown';

r=1:R; Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);

for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    Phir(:,ir)=             (v*diag(besselk(0,d))*v_m1      -v*diag(besseli(0,d))*v_m1*      I0R_m1*K0R)*T_m1*q(1:NLay)';
    Qr(  :,ir)=2*pi*r(ir)*T*(v*diag(besselk(1,d))*v_m1*sqrtA+v*diag(besseli(1,d))*v_m1*sqrtA*I0R_m1*K0R)*T_m1*q(1:NLay)';
end

plotbrug(r,Z,Phir,Qr,s);

%% Solution 720.32 =====================================================

s={
    'Bruggeman (1999) solution 720.32. Steady state.';
    'Fully penetrating wells in all aquifers at r=0.';
    'Constant but different discharges Q. Closed boundary (flux=0) at r=R.'
    };  % phi=phi(x)=drawdown';

r=1:R; Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);

for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    Phir(:,ir)=             (v*diag(besselk(0,d))*v_m1       +v*diag(besseli(0,d))*v_m1*      I1R_m1*K1R)*T_m1*q(1:NLay)';
    Qr(  :,ir)=2*pi*r(ir)*T*(v*diag(besselk(1,d))*v_m1*sqrtA-v*diag(besseli(1,d))*v_m1*sqrtA*I1R_m1*K1R)*T_m1*q(1:NLay)';
end

plotbrug(r,Z,Phir,Qr,s);

%% Solution 720.33 =====================================================

s={
    'Bruggeman (1999) solution 720.33. Steady state.';
    'Fully penetrating wells in all aquifers at r=0 with constant but different discharges Q.';
    'Bunndaries with entrance resistance w to open water at r=R.'
    };  % phi=phi(x)=drawdown';

r=1:R; Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);

F_m1=(sqrtA*I1R+H*I0R)^(-1);
G   =(sqrtA*K1R-H*K0R);
for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    Phir(:,ir)=             (v*diag(besselk(0,d))*v_m1      +v*diag(besseli(0,d))*v_m1*      F_m1*G)*T_m1*q(1:NLay)';
    Qr(  :,ir)=2*pi*r(ir)*T*(v*diag(besselk(1,d))*v_m1*sqrtA+v*diag(besseli(1,d))*v_m1*sqrtA*F_m1*G)*T_m1*q(1:NLay)';
end

plotbrug(r,Z,Phir,Qr,s);

%% Solution 720.34 =====================================================

s={
    'Bruggeman (1999) solution 720.34. Steady state.';
    'All aquifers with fully penetrating cylindrical well screens wih radius R';
    'and constant but different discharges q. Zero flux at r=0.'
    };  % phi=phi(x)=drawdown';

r=1:2*R; Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);

for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    if r(ir)<=R
        Phir(:,ir)=              R*K0R*v*diag(besseli(0,d))*v_m1      *T_m1*q(1:NLay)';
        Qr(  :,ir)=-2*pi*r(ir)*T*R*K0R*v*diag(besseli(1,d))*v_m1*sqrtA*T_m1*q(1:NLay)';
    else
        Phir(:,ir)=             R*I0R*v*diag(besselk(0,d))*v_m1*      T_m1*q(1:NLay)';
        Qr(  :,ir)=2*pi*r(ir)*T*R*I0R*v*diag(besselk(1,d))*v_m1*sqrtA*T_m1*q(1:NLay)';
    end
end

plotbrug(r,Z,Phir,Qr,s);

%% Solution 720.35 =====================================================

% This one needs careful check ! May be Bruggeman constains error.

s={
    'Bruggeman (1999) solution 720.35. Steady state.';
    'Circular polder with radius R, surrounded by a polder wiht different leel, or isolated';
    'reservoir in open water.'
    };  % phi=phi(x)=drawdown';

r=unique([R logspace(0,log10(300),41)]); Nr=length(r);

Phir=zeros(NLay,Nr);
Qr  =zeros(NLay,Nr);

h1=h(1:NLay); h1(2:end)=0;
h2=h(1:NLay); h2(1:end)=0;

%F_m1=(sqrtA*I1R+H*I0R)^(-1);
%G   =(sqrtA*K1R-H*K0R);
for ir=1:length(r)
    [v,d]=eig(r(ir)*sqrtA); v_m1=v^(-1); d=diag(d);
    if r(ir)<=R
        Phir(:,ir)=          (I-R*sqrtA*K1R*v*diag(besseli(0,d))*v_m1)*      h1';
        Qr(  :,ir)=2*pi*r(ir)*T*R*sqrtA*K1R*v*diag(besseli(1,d))*v_m1 *sqrtA*h1';
    else
        Phir(:,ir)=             R*sqrtA*I1R*v*diag(besselk(0,d))*v_m1*      h1';
        Qr(  :,ir)=2*pi*r(ir)*T*R*sqrtA*I1R*v*diag(besselk(1,d))*v_m1*sqrtA*h1';
    end
end

plotbrug(r,Z,Phir,Qr,s);

%%  =====================================================


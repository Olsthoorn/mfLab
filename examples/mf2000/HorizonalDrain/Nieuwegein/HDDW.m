% The problem is tackeled by modeling a 30 m long drain at 5 m below a
% confining layer with no leakage, but with fixed head boundaries at 1000 m
% distance. It is believed that this will yield a proper result for the
% drawdown as in the semi-confined case, the groundwater in the area
% surroundng the drain behave vertiually identical in the case of a semi
% and completely confined situations under the condition that the fixed
% boundary distance applied in the model agrees with the the spreading
% length of the semi-confined system. The here chosen spreading length is
% quite representative for the Dutch situation. The answer is not very
% sensitive to the choice of this distance. The advantage is that we may
% then also apply an analytical solution (see Bruggeman, 1999, p364m
% solution 522.04).
% drain

% This particluar m-file implements Bruggeman's analytical solution for the
% steady state situation. It is a fairly simple one, which, with proper
% appilcation of mirror drains, symmetry and superposition, provides the
% correct solution for the drain under consideration. Exactly the same
% problem can readily be obtained by MODFLOW through mfLab. That sollution
% is implemented in the file mf_adapt in this directory.
%
% Both solution use the same 3D grid to evaluate the drawdown, which
% facilitates comparison.

% The analytical solution is valid for a wll penetrating the roof of a
% confined aqufier of infinite depth. Clearly this is equivalent to our
% (straight) horizontal well in an infinite aquifer with the aquifer roof
% as the mirror plane. If we then turn the well into the horizontal
% position to become a drain of finite extent at a desired distance below
% the roof of the real aquifer, we may implement the roof and the floor of
% the aquifer by application of mirror drains, of which we need an infinite
% number.
%
% Because extraction form a confined aquifer of limited depth has no
% solution, we have to compare the thus computed drawdown with one in a
% reference point, for which we choose a point at the aquifer ceiling at
% distance lambda, for which we took 1000m. The numerical model uses this
% distance for the fixed-head bounary and extent of the model.
%
% TO 100611

%% HDDW specific data

L=30/2;       % [m]   half length of the HDDW
H=50;         % [m]   depth of the aquifer (thickness)
k=15;         % [m/d] hydraulic conductivity
lambda=1000;  %  [m]  distance to reference point
d=1.5;        %  [m]  depth of heart of HDDW below ceiling of aquifer
Q=-19*24/2; rw=0.055;  % [m3/d] Total discharge of (half) the HDDW
xRef=lambda;

close all

%% Grid of points where the drawdown will be computed. The same as of the
% numerical model

xGr= [0:100:1000 0:50:250 0:25:100 0:10:50 0:5:25 0:2:10 0:1:5 0:0.5:2 0.25:1 0:0.125:0.5];
yGr= [0:100:1000 0:50:250 0:25:100 0:10:50 0:5:25 0:2:10 0:1:5 0:0.5:2 0.25:1 0:0.125:0.5];
yGr= [yGr L+(-2:0.5:2) L+(-1:0.25:1)];
zGr=-[0:5:50 0:5:25 0:2:12 0:1:7 d+(-2:25.5:2.25) d+(-0.785:0.25:0.875)];

[xGr,yGr,zGr,xm,ym,zm,DX,DY,DZ,Nx,Ny,Nz]=modelsize3(xGr,yGr,zGr); Z=zGr;


[Xm,Ym,Zm]=meshgrid(xm,ym,zm);  % Generate full 3D coordinate arrays

    s =zeros(size(Xm));        % s will be our drawdown

    for i=-150:150             % 300 mirrors
        
        ZD1=2*(i-1)*H-d; R1=sqrt((Zm-ZD1).^2+Xm.^2);
        ZD2=2*(i-1)*H+d; R2=sqrt((Zm-ZD2).^2+Xm.^2);

        r1Ref=sqrt(ZD1.^2+xRef.^2);
        r2Ref=sqrt(ZD2.^2+xRef.^2);

        ds =+(asinh((Ym+L)./R1  ) -asinh((Ym-L)./R1  ))...
            +(asinh((Ym+L)./R2  ) -asinh((Ym-L)./R2  ))...
            -(asinh((Ym+L)/r1Ref) -asinh((Ym-L)/r1Ref))...
            -(asinh((Ym+L)/r2Ref) -asinh((Ym-L)/r2Ref));
        s=s+ds; err=max(abs(ds(:)));
        fprintf('iter = %4d, ds=%12g s=%12g\n',i,err,mean(s(:)));
        if err<1e-6; break; end
    end
    
    if err>1e-6,
        fprintf('no convergenc in %d iterations\n',i);
    end
    
    s=s*Q/(4*pi*k*L);
    
    fprintf('s=%12g m\n',mean(s(:)));
    
    ccont=-[0:0.02:1];  % contour levels

    Iz=9; % HDDW is in model layer 9 (the analytical solution has no layers)
    
%% Plotting results

subplot(4,1,1); contourf(xm,ym,permute(s(  :,:,Iz),[1,2,3]),ccont); title('xy plane'); xlabel('x'); ylabel('y');
set(gca,'xlim',[0 20],'ylim',[0 20]);

subplot(4,1,2); contourf(xm,zm ,permute(s(end,:,:),[3,2,1]),ccont); title('zx plane'); xlabel('x'); ylabel('z');
set(gca,'xlim',[0 20],'ylim',[-20 0]);

subplot(4,1,3); contourf(ym,zm',permute(s(  :,1,:),[3,1,2]),ccont); title('yz plane'); xlabel('y'); ylabel('z');
set(gca,'xlim',[0 20],'ylim',[-20 0]);

subplot(4,1,4);
PhiWell=s(Iy,Ix,Iz);
plot([-ym(Iy(end:-1:1)) ym(Iy)],[PhiWell(end:-1:1) PhiWell],'b','linewidth',3);
title('Computed head in ideal 30 m long HDDW for 19 m3/h extraction');
xlabel('coordinate along HDDW');
ylabel('drawdown [m]');

mean(PhiWell)

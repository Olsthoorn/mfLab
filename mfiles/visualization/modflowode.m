function h=modflowode(varargin)
% simulates MODFLOW type FD method using ODE solver
% Frans Schaars 151204

[n,~] = getNext(vararging,'double',30);

  %% square grid of nxn , two layers
   gr  =gridObj(0:n,0:n,[1 -3 -7]);
   
   % these paramters seem to be unused, but they are made fields to a struct
   % and passed as such to the function
   p.Kx = gr.const(10);
   p.Ky = gr.const(10);
   p.Kz = gr.const(0.01);
   
   p.CH= gr.const(NaN);  p.CH(:,[1 end],:)=100;
   
   % well in center
   p.Q= gr.const(0);    p.Q( round(n/2), round(n/2), 1)=100;
   
   p.ghc = gr.const(0); p.ghc(9,7)=1000;
   
   p.ghs = gr.const(0); p.ghs(9,7)=90;
   
   p.drc0= gr.const(0); p.drc0(9,7)=0;
   
   p.dre=gr.cont(0);    p.dre(9,7)=0;
   
   p.ric0=gr.const(0);  p.ric0(9,7)=0;
   
   p.rib=gr.const(0.7);
   
   p.ris= gr.Dx;
   
   p.etm=0.00;
   p.ets=gr.const(1);
   p.etx=gr.const(1);
   
   p.S=gr.const(0.5);
   p.S(:,:,2)=0.001;
   
   t=0:.1:10;
   
   p.PHI0=gr.const(100);
   
tspan=t;

%N=rand(1,26)*0.01;
%N=zeros(1,26);
%N(1:end)=0.00;


y0               = p.PHI0(:);
y0(~isnan(p.CH)) = p.CH(~isnan(p.CH));

%opt=odeset('Abstol',1e-3)
tic
[T,Y] = ode15s(@f,tspan,y0,[],p);
toc

h=reshape(Y',[gr.size numel(T)]);

figure

p=surf(h(:,:,1,1)); colorbar

for it=1:numel(t)
    set(p,'zdata' ,h(:,:,1,it));% top layer only
    title(sprintf('heads at time %g d',num2str(T(it))));
    pause(.1);  % wait 0.1 seconds to better see output grahics

    zlim( minmax(h(:)));
    caxis(minmax(h(:))); % adapt caxis
end

%% visualization

V = Y;
V(:,1:end/2)     = V(:,1:end/2).*S(1,1,1);
V(:,end/2+1:end) = V(:,1:end/2).*S(1,1,2);

if nargout==0
    figure;
    subplot(2,2,1);
    plot(T,V);

    subplot(2,2,2);
    plot(T(2:end),diff(sum(V,2)));

    subplot(2,2,3);
    plot(T,Y);

    subplot(2,2,4);
    plot(T(2:end),diff(sum(Y,2)));
end

function dydt = f(t,y,p,gr)
%DYDT -- computes the time derivative of the FD cell heads

    CH =p.CH;
    Q  =p.Q;

    if t>3
        Q=-Q;
    end
    PHI=reshape(y,gr.size);

    %Q(:,:,1) = Q(:,:,1) + interp1(p.t,p.N,t,'linear','extrap').*Dx(:,:,1).*Dy(:,:,1);

    %conductance in x,y, en z-richting
    Rx = 0.5* gr.DX ./(gr.DY .* gr.DZ) ./ p.Kx;
    Ry = 0.5* gr.DY ./(gr.DX .* gr.DZ) ./ p.Ky;
    Rz = 0.5* gr.DZ ./(gr.DX .* gr.DY) ./ p.Kz;
    
    Cx = 1 ./ (Rx(:,1:end-1,:) + Rx(:,2:end,:));    
    Cy = 1 ./ (Ry(1:end-1,:,:) + Ry(2:end,:,:));
    Cz = 1 ./ (Rz(:,:,1:end-1) + Rz(:,:,2:end));

    Qx(:,2:end-1,:) = -diff(PHI,1,1).*Cx;
    Qy(2:end-1,:,:) = -diff(PHI,1,2).*Cy;
    Qz(:,:,2:end-1) = -diff(PHI,1,3).*Cz;

    dydt=zeros(gr.size);

    dydt(1:end-1,:,:)=dydt(1:end-1,:,:)-Qx;
    dydt(:,1:end-1,:)=dydt(:,1:end-1,:)-Qy;
    dydt(:,:,1:end-1)=dydt(:,:,1:end-1)-Qz;

    dydt(2:end,:,:)=dydt(2:end,:,:)+ Qx;
    dydt(:,2:end,:)=dydt(:,2:end,:)+ Qy;
    dydt(:,:,2:end)=dydt(:,:,2:end)+ Qz;

    dydt=(Q(:)+dydt(:))./p.S(:);

    dydt(~isnan(CH))=0;
    fprintf(' %g',t);
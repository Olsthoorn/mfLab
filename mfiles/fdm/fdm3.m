function [Phi,Q,Qx,Qy,Qz,Psi]=fdm3(gr,kx,ky,kz,IBOUND,IH,FQ) %% future: ,hdr,Cdr,hrv,hbr,Crv,hgh,Cgh)
%FDM3 a 3D block-centred steady-state finite difference  with selftest
%
% Example:
%    [Phi,Q,Qx,Qy,Qz,Psi]=fdm3(x,y,Z,x,kx,ky,kz,IBOUND,IH,FQ)
%    fdm3() % runs selftest in axial symmetric mode
%
% INPUT:
%    x,y,Z mesh coordinates
%    Z may be a vector, a 3D vertical vector or a full size 3D array of size Ny,Nx,Nz+1
%    kx,ky,kz conductivity arrays
%    IBOUND as in MODFLOW
%    IH=initial heads (NaN for ordinary points),
%    FQ=fixed nodal flows
% OUTPUT:
%    Phi,Q computed heads and cell balances
%    Qx,Qy,Qz is cell face flow, positive along positive axis direction
%    Psi   stream function, only in axial symmetric mode
%    Future:  hdr=elevation drains, Cdr=resistance drains
%             hrv=elevation river , Crv=resistance river, hbr=elevation bottom river
%             hgh=elevation general head, Cgh is resistance general head
%
% See also: fmd2t fdm2c fdm2ct fdm3 fdm3t
%
% TO 991017  TO 000530 001026 070414 070513
% TO 080226 implemented inactive cells and true fixed heads
% TO 090216 small edits

% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin<7
    selftest()
    return
end

Nodes = reshape(1:gr.Nxyz,gr.size);               % Node numbering
IW=Nodes(:,1:end-1,:);    IE=Nodes(:,2:end,:);
IN=Nodes(1:end-1,:,:);    IS=Nodes(2:end,:,:);
IT=Nodes(:,:,1:end-1);    IB=Nodes(:,:,2:end);

% resistances and conducctances
warning('off','all');  % to allow zero conductivities for inactive cells
if gr.AXIAL
    fprintf('Running %s in axial symmetric mode\n',mfilename);
    RX = bsxfun(@times, log(gr.xm( :,2:end  )./gr.xGr(2:end-1)), 1./ (2 * pi * kx(:,2:end  ,:) .* gr.DZ(:,2:end  ,:))) + ...
         bsxfun(@times, log(gr.xGr(:,2:end-1)./gr.xm( 1:end-1)), 1./ (2 * pi * kx(:,1:end-1,:) .* gr.DZ(:,1:end-1,:)));
    RZ = 0.5 * gr.DZ ./ bsxfun(@times, pi * (gr.xGr(2:end).^2 - gr.xGr(1:end-1).^2), kz);
    Cx = 1./RX;
    Cy = zeros(gr.Ny-1,gr.Nx,gr.Nz);
    Cz = 1./(RZ(:,:,1:end-1) + RZ(:,:,2:end));    
else
    fprintf('Running %s in linear mode\n',mfilename);
    RX=0.5*gr.DX./(gr.DY.*gr.DZ.*kx); Cx=1./(RX(:,1:end-1,:)+RX(:,2:end,:));
    RY=0.5*gr.DY./(gr.DZ.*gr.DX.*ky); Cy=1./(RY(1:end-1,:,:)+RY(2:end,:,:));
    RZ=0.5*gr.DZ./(gr.DX.*gr.DY.*kz); Cz=1./(RZ(:,:,1:end-1)+RZ(:,:,2:end));
end
warning('on','all');

A=sparse([IE(:);IW(:);IN(:);IS(:);IT(:);IB(:)],...
         [IW(:);IE(:);IS(:);IN(:);IB(:);IT(:)],...
        -[Cx(:);Cx(:);Cy(:);Cy(:);Cz(:);Cz(:)],...
         gr.Nxyz,gr.Nxyz,7*gr.Nxyz);                 % System matrix
Adiag= -sum(A,2);                               % Main diagonal

IAct = IBOUND~=0; % active cells
I    = IBOUND >0; % active cells but not fixed heads = cells with heads to be computed
Ifh  = IBOUND <0; % active cells with fixed heads


Phi=IH(:); FQ =FQ(:); Q  =FQ;

if any(I(:))  % something to compute
    Phi( I)=spdiags(Adiag(I   ),0,A(I   ,I   ))\(FQ(I)-A(I,Ifh)*Phi(Ifh)); % solve
end
if any(IAct)
    Q(IAct)=spdiags(Adiag(IAct),0,A(IAct,IAct))*Phi(IAct );		           % reshape
end

Phi= reshape( Phi,size(IBOUND));
Q  = reshape( Q  ,size(IBOUND));

if nargout>2
    Qx=-Cx.*diff(Phi,1,2)*sign(gr.xGr(end)-gr.xGr(1)); Qx(isnan(Qx))=0;  % Flow across horizontal cell faces
    Qy=-Cy.*diff(Phi,1,1)*sign(gr.yGr(end)-gr.yGr(1)); Qy(isnan(Qy))=0;  % Flow across vertical cell faces
    Qz=-Cz.*diff(Phi,1,3)*sign(gr.Z(end)-gr.Z(1));     Qz(isnan(Qz))=0;  % Flow across vertical cell faces
    
    if gr.AXIAL
        Psi = zeros(gr.Ny,gr.Nx-1,gr.Nz+1);
        Psi(:,:,1:end-1) = Qx;
        Psi = cumsum(Psi(:,:,end:-1:1),3);
        Psi = Psi(:,:,end:-1:1);
    else
        Psi = NaN;
    end
end

function selftest()
    %% selftest() --- tests fdm3 in axial symmetric mode
    
    fprintf('selftest %s for axial symmetric mode\m',mfilename);
    
    % Grid
    xGr = logspace(0,4,51);
    yGr = [-0.5 0.5];
    zGr = 0:-1:-20;
    gr = gridObj(xGr,yGr,zGr,'AXIAL',true); % axial mode
    
    % Fixed head locations (top)
    IBOUND = gr.const(1);
    IBOUND(:,:,1)=-1;
    IBOUND(:,end,:)=-1;
    
    % Fixed flows
    Qw     = -2400;
    zW(2)  = gr.zGr(end)+0.2*diff(gr.zGr([end 1])); % 20% of elevation
    zW(1)  = gr.zGr(end)+0.6*diff(gr.zGr([end 1])); % 80% of elevation
    Iw     = gr.Nx * gr.Ny * (find(gr.zm>zW(2) & gr.zm<zW(1)) - 1) + 1;     % well in first row first column only
    FQ     = gr.const(0);
    FQ(Iw) = Qw * gr.DZ(Iw)/sum(gr.DZ(Iw));
    
    % Initial heads
    IH     = gr.const(0);
    
    % Conductivities
    kx     = gr.const(10);
    %top layer resistance 1000 d
    kx(:,:,1) = gr.DX(:,:,1)/1000/2;
    
    ky     = kx;
    kz     = kx;
    
    kz(Iw) = 1000; % high vertical conductivity in well screen

    % 3D model (axial symmetric mode
    [Phi,Q,Qx,Qy,Qz,Psi]=fdm3(gr,kx,ky,kz,IBOUND,IH,FQ);

    % Contouring
    hrange = ContourRange(Phi,50); % head
    prange = ContourRange(Psi,50); % stream function

    figure('pos',screenPos(0.6));
    ax = axes('nextPlot','add','clim',hrange([1 end]),'xscale','log');

    title(ax,['self test ' mfilename]); xlabel('r [m]'); ylabel('z [m]');
    contourf(ax,gr.xm,gr.zc,XS(Phi),hrange);
    contour(ax,gr.xp,gr.zp,XS(Psi),prange,'color',grey)
    hb = colorbar();
    hb.Label.String='head [m]';

    % Print overview of outcomes
    fprintf('Phi    %12g ... %12g\n',min(Phi(:)),max(Phi(:)));
    fprintf('Q      %12g ... %12g, sum Q: in = %12g sum out = %12g\n',min(Q(:)),max(Q(:)),sum(Q(Q>0)),sum(Q(Q<0)));
    fprintf('Qx     %12g ... %12g\n',min(Qx(:)),max(Qx(:)));
    if ~gr.AXIAL
        fprintf('Qy     %12g ... %12g\n',min(Qy(:)),max(Qy(:)));
        set(ax,'xscale','log')
    end
    fprintf('Qz     %12g ... %12g\n',min(Qz(:)),max(Qz(:)));
    fprintf('Psi    %12g ... %12g\n',min(Psi(:)),max(Phi(:)));

    fprintf('done !\n');

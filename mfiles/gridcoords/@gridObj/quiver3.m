function h = quiver3(o,B,plane,I,S)
% h = gridObj.quiver3(B,plane,I) -- plot specific discharge vectors implied by
% B(end), plane and I
% B originates from readBud. Always use last. Use B(it) in call if it is required
% plane is 'zx' I=iy, 'zy' I=ix, or 'xy' I=iz
% TO 120827

if nargin<5, S=0; end  % no scaling of the errors

Ix = 1:o.Nx;
Iy = 1:o.Ny;
Iz = 1:o.Nz;

Qx = B(end).term{strmatchi('FLOWRIGHTFACE',B(end).label)};
Qy = B(end).term{strmatchi('FLOWFRONTFACE',B(end).label)};
Qz = B(end).term{strmatchi('FLOWLOWERFACE',B(end).label)};

I=I(1); 

switch plane
    case {'zx','x'}
        Iy=I;
    case {'zy','y'}
        Ix=I;
    case 'xy'
        Iz=I;
end

Qx = cat(2,Qx(:,1,:),Qx);
qx = ((Qx(:,1:end-1,:)+Qx(:,2:end,:))/2)./o.DZ./o.DY;

Qy = cat(1,Qy(1,:,:),Qy);
qy = ((Qy(1:end-1,:,:)+Qy(2:end,:,:))/2)./o.DZ./o.DX;

Qz = cat(3,Qz(:,:,1),Qz);
qz = ((Qz(:,:,1:end-1)+Qz(:,:,2:end))/2)./o.DX./o.DY;

S=1;

h = quiver3(o.XM(Iy,Ix,Iz),o.YM(Iy,Ix,Iz),o.ZM(Iy,Ix,Iz),qx(Iy,Ix,Iz),qy(Iy,Ix,Iz),qz(Iy,Ix,Iz),S);


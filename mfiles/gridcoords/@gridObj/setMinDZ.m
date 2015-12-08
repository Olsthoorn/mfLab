function [HK,VK,o] = setMinDZ(o,HK,VK,MINDZ,layers)
% gr.setMinDZ(MINDZ,layers)
% Simple minded function to set the minimum thickness of layers to MINDZ.
% Thickness is simply asserted at the bottom of the layer, wihile
% disregarging the rest of the model.
%
% MINDZ is important to reduce computation time of SEAWAT and MT3DMS
% 
% TO 120612

if any(o.LAYCBD)~=0
    error(['%s: Sorry, byt this procdure is only valid if all LAYCBD==0,\n',...
          'first remove your CBD layers with gridOb.removeCBD'],mfilename);
end

if nargin<4, MINDZ = o.MINDZ; end
if nargin<5, layers = 1:o.Nlay; end

errstr = '%s size %s=[%d %d %d] must match size of grid=[%d %d %d]';

if ~all(size(HK)==o.size), error(errstr,mfilename,'HK',size(HK),gr.size); end
if ~all(size(VK)==o.size), error(errstr,mfilename,'VK',size(VK),gr.size); end

if MINDZ<0, error('%s: MINDZ=%g must be greater than zero',mfilename,MINDZ); end

layers = min(o.Nlay-1,max(1,layers));

N  = o.Nx*o.Ny;
Z  = o.Z;

m = NaN(o.Nlay,length(layers));

for iz=layers
    m(:,iz) = XS(min(min(-diff(Z,1,3))))
    DZ = Z(:,:,iz)-Z(:,:,iz+1);

    I = find(DZ<MINDZ);
    if ~isempty(I)
        Itop = I+(iz-1)*N;

        % Correction necessary or the plane below the layer bottom Z(:,:,Itop+2*N)
        dz   = MINDZ - DZ(I);

        % Make sure that the new HK for this part of the layers is unaltered
        HK(Itop)  = (HK(Itop).* DZ(I) + HK(Itop+N).*dz)/MINDZ;

        % Make sure that the new VK for this part of the layer is unchanged
        VK(Itop)  = MINDZ./(DZ(I)./ VK(Itop) + dz./VK(Itop+N));

        % The correction of the underlying layer is implied by the change of
        % the vertical position of the bottom of the current layer,  i.e. the
        % top of the underlying layer
        Z (Itop+  N) = Z(Itop  )-MINDZ;
        for i=iz+1:o.Nz
            Z (I+i*N) = min(Z(I+i*N),Z(I+(i-1)*N)-o.MINDZ);
        end
    end
end

o = gridObj(o.xGr,o.yGr,Z,o.LAYCBD);

% if islogical(layers), layers=find(layers); end
% 
% layers = layers(layers>=1 & layers<= o.Nz);
% 
% DZlayers = o.DZ(:,:,layers);
% 
% DZlayers(DZlayers<MINDZ)=MINDZ;
% 
% DZ=o.DZ;
% 
% DZ(:,:,layers)=DZlayers;
% 
% Z=o.Z;
% 
% for iz=1:size(DZ,3);
%     Z(:,:,iz+1)=Z(:,:,iz)-DZ(:,:,iz);
% end
% 
% gr = gridObj(o.xGr,o.yGr,Z,o.LAYCBD,o.MINDZ,o.AXIAL);

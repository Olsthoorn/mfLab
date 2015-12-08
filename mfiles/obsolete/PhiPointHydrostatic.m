function phip=PhiPointHydrostatic(Phi,Z,Rho)
% PhHIPOINTHYDROSTATIC hydrostatic point-water head at centers given by Z
%    Phi,Z and Rho are all 3D arrays,
%    Z may have one row more than Rho it then is the elevation of the tops
%    and bottoms of all layers
%    Rho may be absolute or relative Rho=(rho-rhof)/rhof 
%    TO 090308
%
%    EXAMPLE
%       phip=PHIPOINTHYDROSTATIC(Phi,z,rho)
%
% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

[Ny.Nx,Nz]=size(Phi);
if size(Z,3)>Ny,
    Z = 0.5*(Z(:,:,1:end-1)+Z(:,:,2:end));
end

P1  = repmat( (Phi(:,:,1)-Z(:,:,1)) .* Rho(:,:,1),[1,1,Nz]);  % gravity is immaterial
Pm  = P1-cumtrapz(Z,Rho,3);

phip = Pm./rho+Z;

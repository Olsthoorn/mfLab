function [dx,dy,r] = dist(gp1,gp2)
% [dx,dy,r] = dist(gp1,gp2);  computes distance in m between two googleMap points
% [r      ] = dist(gp1,gp2);  computes only the radial istance on the globe
%
% TO 110501 121008

R2 = 6371007.2;  % The Earth's authentic radius (Geodetic Union);

gp1=gp1(1);

if nargin<2
    [Lon2 Lat2]=ginput;
    gp2 = googlePointObj(Lat2,Lon2);
end

dLam = [gp2.lam]-gp1.lam;
dPhi = [gp2.phi]-gp1.phi;

a = [cos( gp1.phi)  *cos( gp1.lam)  ...
     cos( gp1.phi)  *sin( gp1.phi)  ...
     sin( gp1.phi) ];
b = [cos([gp2.phi]).*cos([gp2.lam]) ...
     cos([gp2.phi]).*sin([gp2.lam]) ...
     sin([gp2.phi])];

 for i=numel(gp2):-1:1
     r(i) = R2 * acos(dot(b(i,:),a));
 end
 
if nargout>1
     dx = R2 * cos(([gp2.phi] + gp1.phi)/2) .* dLam;
     dy = R2 * dPhi;
end

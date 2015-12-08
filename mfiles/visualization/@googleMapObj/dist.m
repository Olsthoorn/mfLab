function [dx,dy,r] = dist(gp1,gp2)
% [dx,dy,r] = dist(gp1,gp2) compute distance in m between two googleMap points
% [dx,dy,r] = gp1.dist(gp2) compute distance in m between two googleMap points
% [r      ] = gp1.dist(gp2) compute radial distance if nargout = 1
% gp1 and gp2 are googlePointObj's
%
% TO 110501 121008

R2 = 6371007.2;  % The Earth's authentic radius (Geodetic Union);

if nargin<2
    [Lon Lat]=ginput;
    gp2 = googlePointObj(Lat,Lon);
end

dlam = [gp2.lam] - gp1.lam;
dphi = [gp2.phi] - gp1.phi;

a = [cos(gp1.lam)*cos(gp1.phi) sin(gp1.lam)* cos(gp1.phi)  sin(gp1.phi)];
b = [cos(lam(:)).cos(phi(:)) sin(lam(:)).*cos(phi(:)) sin(phi(:))];

for i=numel(gp2):-1:1
    dx(i) = R2 * dlam * cos((gp2(i).phi + gp1.phi)/2);
    dy(i) = R2 * dphi;
    r(i)  = R2 * acos(dot(a,b(i,:)));
end



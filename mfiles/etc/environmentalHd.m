function hEnv = environmentalHd(z,hd,salinity)
%ENVIRONMENTALHD computes environmental head (Seawat boundary conditions)
%
% USAGE:
%    hEnv = environmentalHd(z,hd,salinity) -- computes environmental head
%
% Assuming salinity of ocean water = 1
%
% z          is elevation of model cells of one column
% hd         is environmental head in that column
% salinity   is relative concentration (fresh=0 and saline=1), so that
% salinity is proportional to Rhos-Rhof
%
%    All cells in column are assumed to be active
%
%    assume z=3D, hd=3D or scalar, salinity=3D rho/rhof = 1+ddeltadc*salinity
%    where salinity is relative to ocean water hence ocean water  = 1.
%
% used as an alternative to or for verificatin of CHDDENSOPT for CHD when
% used with Seawat to specify a hydrotatic vertical boundary of the model.
% See Langevin et al. (2006) manual for SEAWAT version 4.
%
% TO 120510

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

z        =  XS(z);
salinity =  XS(salinity);
ddeltadc = 0.025;
rho = 1+ ddeltadc * salinity;

zt = z(1:end-1,:);
zb = z(2:end,  :);
zm = (zt+zb)/2;
dz = abs(diff(z,1,1));

hd =  ones(size(zt(:,1))) * hd;


%% See if head is above top of model
Dz0rho =  ones(size(zt(:,1)))* ((hd(1,:)-zt(1,:)) .* (hd(1,:)>zt(1,:)) .* rho(1,:));

Dz = dz.*(hd>zt) + (hd-zb).*(hd<=zt & hd>=zb);

SumRhoDz = Dz0rho + cumsum(Dz.*rho) -rho.*dz/2;

hEnv = zm + SumRhoDz./rho;

hEnv = XS(hEnv);

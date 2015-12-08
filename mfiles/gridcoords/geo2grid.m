function kGr=geo2grid(zGr,zGeo,kGeo)
%GEO2GRID maps geology given by zGeo and kGeo to grid given by zGr
%
% Example:
%   kGr=geo2grid(zGr,zGeo,kGeo)
%
%   zGeo are layer interfaces Ngeo+1
%   kGeo are layer conductivities Ngeo
%   zGr  are model layer elevations Nz+1
%
% See also: gridsTransfer
%
%   TO 091216

zGr = zGr(:); zGeo=zGeo(:); kGeo=kGeo(:);
zm  = 0.5*(zGr(1:end-1)+zGr(2:end));
dz  = abs(diff(zGr));

% merge the two grids
Z   = unique([zGr(:); zGeo(:)]); Z=Z(end:-1:1);
Zm  = 0.5*(Z(1:end-1)+Z(2:end));
dZm = abs(diff(Z));

% find which geo layer matches the merged grid
I=zeros(size(Zm));  % I gets the geo layer numbers
for i=1:length(zGeo)
    I=I+(Zm<zGeo(i));
end
kDmerge=kGeo(I).*dZm;  % merged layers are always unique

% find which mesh layer matches the merged layer
J=zeros(size(Zm));
for j=1:length(zm)
    J=J+(Zm<zGr(j));
end

% grab together the matched merged layers
kDGr=NaN(size(zm));
for i=1:length(zm)
    kDGr(i)=sum(kDmerge(J==i));
end
kGr=kDGr./dz;

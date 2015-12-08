function [Psi]=streamf(xMesh,yMesh,QRF,QFF)
% [Psi]=streamf(xMesh,yMesh,QRF,QFF)
% compute stream function using QRF, QFF is not used
% Psi contains values per node. All water extraction is assumed to flow through top of system.
% TO 010811

[Ny,Nx]=size(QRF);
if Nx<length(xMesh),
   QRF=[zeros(Ny,1),QRF,zeros(Ny,1)];
end

[Ny,Nx]=size(QFF);
if Ny<length(yMesh);
   QFF=[zeros(1,Nx);QFF;zeros(1,Nx)];
end


Psi=flipud(cumsum(flipud([QRF;zeros(size(xMesh))])));

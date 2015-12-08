function [Phi,Q,Psi,Qx,Qy]=fdm2(xGr,yGr,Kx,Ky,IBOUND,IH,FQ,varargin)
%FDM2 a 2D block-centred steady-state finite difference model
%
% Example:
%    [Phi,Q,Psi,Qx,Qy]=fdm2d(xGr,yGr,Kx,[Ky],IBOUND,FH,FQ [,radial])
%
% Inputs:
%    xGr,yGr mesh coordinates, Kx,Ky conductivities, if Ky==[], Ky=Kx
%    xGr must be a hor  vector or a grid
%    yGr must be a vert vector or a grid
%    IBOUND(Ny,Nx): Inactive cells IBOUND==0, fixed heads IBOUND<0
%    IH=initial heads
% Outputs:
%    Phi,Q computed heads and cell balances
%    Qx(Ny-1,Nx-2) horizontal cell face flow positive in positive xGr direction
%    Qy(Ny-2,Nx-1) vertial    cell face flow, postive in positive yGr direction
%    Psi(Ny,Nx-2)  stream function
%
% See also: fmd2t fdm2c fdm2ct fdm3 fdm3t
%
% TO 991017  TO 000530 001026 070414 090314 101130

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if isempty(varargin)
    [Phi,Q,Psi,Qx,Qy]=fdm2c(xGr,yGr,Kx,[],Ky,IBOUND,IH,FQ);
else
    [Phi,Q,Psi,Qx,Qy]=fdm2c(xGr,yGr,Kx,[],Ky,IBOUND,IH,FQ,varargin{1});
end
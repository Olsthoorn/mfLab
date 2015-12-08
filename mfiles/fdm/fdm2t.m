function [Phi,Qt,Qx,Qy,Qs]=fdm2t(x,y,t,kx,ky,S,IBOUND,IH,FQ,varargin)
%FDM2T a 2D block-centred transient finite difference model
%
% Example:
%    [Phi,Q,Qx,Qy,Qs]=fdm2dt(x,y,t,kx,ky,S,IBOUND,IH,FQ [,radial])
%
% INPUT:
%  x(Nx+1)         x-coordinate of mesh/grid
%  y(Ny+1)         y-coordinate of mesh/grid
%  t(Nt+1)         time for output. 0 will be added if necessary
%  kx(Ny,Nx)       conductivity in x-direction
%  ky(Ny,Nx)       same in y direction, ky=kx if ky==[]
%  S(Ny,Nx)        primary storage (S+Sy)
%  IBOUND(Ny,Nx)   same as in MODFLOW;
%  IH(Ny,Nx)       initial head
%  FH(Ny,Nx)       fixed heads (NaN for ordinary points), Q=fixed nodal flows
%  Radial          Arbitrary input caused model to assume axial flow with
%                 x=r and FQ are ring flows while als flows have dimension
%                 L3/T instead of L2/t. For clearness use 'R' or 'Radial'
%                 as input at the function call.
%
% OUTPUT
%  Phi(Ny,Nx,Nt+1) computed heads with Phi(Ny,Nx,1) initial heads for t=0
%, Qt(Ny,Nx,Nt)    computed total cell balance during timestep it
%  Qx(Ny,Nx-1,Nt)  hor.  cell face flow in x-direction positive along increasing col number
%  Qy(Ny-1,Nx,Nt)  vert. cell face flow in y-directin  positive along increasing row number
%  Qs(Ny,Nx,Nt)    storage change of node during timestep it
%
% See also: fmd2t fdm2c fdm2ct fdm3 fdm3t
%
% TO 991017  TO 000530 001026 070414 080301

% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%% Allow simple input

%% call fdm2ct
if isempty(varargin)
    [Phi,Qt,Qx,Qy,Qs]=fdm2ct(x,y,t,kx,[],ky,S,IBOUND,IH,FQ);
else
    [Phi,Qt,Qx,Qy,Qs]=fdm2ct(x,y,t,kx,[],ky,S,IBOUND,IH,FQ,varargin{1});
end


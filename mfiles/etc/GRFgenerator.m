function X=GRFgenerator(x,y,rho,sigma)
%GRFGENERATOR Gaussian Random field generator according to Hoo (2004)
%
% USAGE:
%    X=GRFgenerator(x,y,rho,sigma)
% 
% It is perfectly similar to a finite difference model.
% Refer to paper for further documentation or try it.
%
% Useful for testing geophysics.
%
% TO 070502

% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later


x=x(:)'; Nx=length(x); %dx=diff(x); xm=0.5*(x(1:end-1)+x(2:end));
y=y(:);  Ny=length(y); %dy=abs(diff(y));

Nodes = reshape(1:Nx*Ny,Ny,Nx);               % Node numbering
IE=Nodes(:,2:end);   IW=Nodes(:,1:end-1);
IS=Nodes(2:end,:);   IN=Nodes(1:end-1,:);

% resistances and conducctances
Cx=0.25*ones(Ny,Nx-1); Cy=0.25*ones(Ny-1,Nx);

A=sparse([IE(:);IW(:);IN(:);IS(:)],...
         [IW(:);IE(:);IS(:);IN(:)],...
         [Cx(:);Cx(:);Cy(:);Cy(:)],...
         Ny*Nx,Ny*Nx,5*Ny*Nx);                 % System matrix
Adiag=-sum(A,2);                               % Main diagonal

e=randn(Ny,Nx)*sigma;                          % Normal random field with zero correlation and std sigma

X=reshape(spdiags(Adiag,0,rho*A)\e(:),Ny,Nx);                 % solve using spatial smoothing factor rho

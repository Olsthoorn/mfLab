function [FQ]=flowbound(Q,zMax,zmin,Col,z,FQ)
%[FQ]=flowbound(Q,zMax,zmin,Col,z,FQ) adds flow Q in screen between zMax and zmin to nodes in column Col
% q=Q/length(screen) = Q/(zMax=zmin)
% TO 990419

if z(end)<z(1), z=-z; zMax=-zMax; zmin=-zmin; end

if zmin>zMax; dum=zmin; zmin=zMax; zMax=dum; end

q=Q/(zMax-zmin);

zmid=(z(1:end-1,Col)+z(2:end,Col))./2;

% right hand side vector
FQ(:,Col)= q.*...
   ([max(0,(min(zMax,zmid)-max(zmin,z(1:end-1,Col))));0]+[0;max(0,(min(zMax,z(2:end,Col))-max(zmin,zmid)))]);

function A = sysmat(kD,c)
% system matrix multi-layer solutions

nlay = numel(kD);
kD = kD(:);
c  = c(:);

A=  -diag( 1./(kD(2:nlay  ).*c(2:nlay)),-1)+...
    +diag( 1./(kD(1:nlay  ).*c(1:nlay))+1./(kD(1:nlay).*[c(2:nlay); Inf]), 0)...
    -diag( 1./(kD(1:nlay-1).*c(2:nlay)),1);

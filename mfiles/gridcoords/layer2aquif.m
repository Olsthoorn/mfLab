function iaq=layer2aquif(iz,LAYCBD,Nz)
%LAYER2AQUIF converts layer number (LAY+CBD) to aquifer number (LAY) using full LAYCBD
%
% Example:
%     iaq=layer2aquif(iz,LAYCBD,Nz)
%
% where Nz = layers + cbdlayers
%
% See also: gridObj
%
% TO 120407

IAq =cumsum(isAquifer(Nz,LAYCBD));
iaq=IAq(iz);

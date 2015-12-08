function z = zRel2z(o,zRel,Idx)
%%ZREL2Z -- converts zRel to Z (relative z-coordinates to absolute z
% coordinates)
%
% USAGE z = gridObj.zRel2z(zRel)
%
% zRel is relative z-coordinates which starts at 0 at the top of the grid,
% is 1 at the bottom of the first layer etc and is gr.Nlay at the bottom of
% the system. Doe not work yet for grids with confining beds.
%
% TO 131021

Layer = max(1,ceil(zRel(:)));
frac  = zRel(:) - Layer + 1;
Layer = min(Layer,o.Nlay);

Idx = o.IdTop(Idx(:)) +(Layer-1)*o.Nxy;

z = o.Z(Idx) - frac .* o.DZ(Idx);
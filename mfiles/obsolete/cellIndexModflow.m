function Idx=(well,gr,frac)
%CELLINDEXMODFLOW produces global cell indices for a well screen point for MODFLOW grid
%
% to get the cellIndex for the MODFLOW grid, to be used as recirculation
% wells with MT3DMS or SEAWAT use method wellObj.screenPoint(gr,frac)
%
% wellObj.screenPiont(gr,fr) finds the cell of the well that matches the
% fraction of the cell length computed from its top. It adds values to the
% wells UserData, one of which is idxMT3D the global index in the
% MODFLOW/MT3D grid where the well screen point is. Retrieve like this
%
% well = well.screenPoint(gr,frac);
% well(iw).UserData.screenPoint.idxMT3D;
% see
% well(iw).UserData.screenPoint to view additonal stored data.
%
% SEE ALSO: wellObj.screenPoint(gr,frac)
%
% TO 090317 090509 100521 100830 121127
%
% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

well = well(1).screenPoint(gr,frac);
Idx = well.UserData.screenPoint.idxMT3D;

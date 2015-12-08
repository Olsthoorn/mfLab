function o=setWell(o,gr,HK) 
%% well = well.setWell(gr,HK) places the well in the grid
% by computing indices and fraction of the discharge from the
% cells penetrated by the well screen.
%
% TO 120512

% now replaced by
o = o.toGrid(gr,HK);



function ZETA = zeta(gr,surfaces)
% ZETA = gr.zeta(surfaces) --- obsolete !! (per 140910)
% where surfaces is array [Ny,Nx,Nsurf] with the elevation of the zeta
% surfaces required by the SWI package
%
% Now obsolete, use surface direct in mf_adapt
mindz = 0.001;

for iSurf = size(surfaces,3):-1:1
    ZETA(iSurf).values = gr.const(NaN);
    for iLay = 1:gr.Nlay
        ZETA(iSurf).values(:,:,iLay) = ...
            max(gr.ZBlay(:,:,iLay)+mindz, ...
                min(gr.ZTlay(:,:,iLay)-mindz,surfaces(:,:,iSurf)) ...
                );
    end
end

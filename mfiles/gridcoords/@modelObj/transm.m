function [Tran,HK] = transm(o,Ix,Iy,Iz)
% [Tran,HK] = modelObj.tran([Ix [,Iy [,Iz]]])- compute transmissivity of layers
% for Ix,Iy,Iz.
% TO 120829

iHK = strmatchi('HK',{o.name},'exact');

if ~iHK
    error('%s: HK not in Model');
else
    HK= o(iHK).var;
end

gr = o.grid();

if nargin<4, Iz=1:gr.Nz; end
if nargin<3, Iy=1:gr.Ny; end
if nargin<2, Ix=1:gr.Nx; end

DZ  = gr.DZ(Iy,Ix,Iz);
HK  = HK(Iy,Ix,Iz);
Tran= DZ.*HK;
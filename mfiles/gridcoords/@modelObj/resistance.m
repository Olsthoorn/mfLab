function [c,VK] = resistance(o,Ix,Iy,Iz)
% [c,VK] = modelObj.tran([Ix [,Iy [,Iz]]])- compute resistance of layers
% for Ix,Iy,Iz.
%
% TO 120829

iHK = strmatchi('HK',{o.name},'exact');

if ~iHK
    error('%s: HK not in Model');
else
    VK= o(iHK).var;
end

gr = o.grid();

if nargin<4, Iz=1:gr.Nz; end
if nargin<3, Iy=1:gr.Ny; end
if nargin<2, Ix=1:gr.Nx; end

DZ  = gr.DZ(Iy,Ix,Iz);
VK  = VK(Iy,Ix,Iz);
c= DZ./VK;
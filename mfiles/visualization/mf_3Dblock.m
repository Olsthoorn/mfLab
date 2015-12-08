function h=mf_3Dblock(varargin)
%MF_3DBLOCK plots a 3D block with cutout defined by ix,iy,iz
%
% USAGE:
%   h=mf_3Dblock(XM,YM,ZM,C,ix,iy,iz[,facealpha[,wireframeclr]);
%
% Example:
%    h=mf_3Dblock(XM,YM,ZM,C,ix,iy,iz,facealpha,clr)
%    h=mf_3Dblock(ax,XM,YM,ZM,C,ix,iy,iz,facealpha,clr)
%    h=mf_3Dblock(grid,ix,iy,iz,facealpha,clr)
%    h=mf_3Dblock(ax,grid,ix,iy,iz,facealpha,clr)
%
%   h is a 3 by 6 matrix with handles to all 3x6 block surfaces
%   use -ix, -iy, iz to get the block ix:Nx, iy:Ny, iz:Nz respectively.
%
%   For normal view(3) use ix,-iy,-iz
%
%   use ix>=Nx, iy>=Ny iz>Nz to get total block without cut 
%
% See also: mf_wireframe
%
% TO 110501

hold on; % just to make sure we will see all

h=varargin{1};
if isaxis(varargin{1})
    ax=h;
    varargin(1)=[];
else
    ax = gca;
end

set(ax,'nextplot','add');

if strcmpi(class(h),'gridObj')
    gr = varargin{1};
    varargin(1)=[];
else
    gr=gridObj(varargin{1},varargin{2},varargin{3});
    varargin(1:3)=[];
end

C         = varargin{1}; varargin(1)=[];
ix        = varargin{1}; varargin(1)=[];
iy        = varargin{1}; varargin(1)=[];
iz        = varargin{1}; varargin(1)=[];
faceAlpha = varargin{1}; varargin(1)=[];
clr       = varargin{1};

if ix>0, Ix=1:ix; else Ix=abs(ix):gr.Nx; end; ix=min(abs(ix),gr.Nx);
if iy>0, Iy=1:iy; else Iy=abs(iy):gr.Ny; end; iy=min(abs(iy),gr.Ny);
if iz>0, Iz=1:iz; else Iz=abs(iz):gr.Nz; end; iz=min(abs(iz),gr.Nz);

mf_wireframe(ax,...
             gr.XM([1 iy gr.Ny],[1 ix gr.Nx],[1 iz gr.Nz]),...
             gr.YM([1 iy gr.Ny],[1 ix gr.Nx],[1 iz gr.Nz]),...
             gr.ZM([1 iy gr.Ny],[1 ix gr.Nx],[1 iz gr.Nz]),...
             clr);
         
h=NaN(3,6);

h(1,:)=mybox(ax,gr.XM(:,Ix,:),gr.YM(:,Ix,:),gr.ZM(:,Ix,:),C(:,Ix,:),faceAlpha);
h(2,:)=mybox(ax,gr.XM(Iy,:,:),gr.YM(Iy,:,:),gr.ZM(Iy,:,:),C(Iy,:,:),faceAlpha);
h(3,:)=mybox(ax,gr.XM(:,:,Iz),gr.YM(:,:,Iz),gr.ZM(:,:,Iz),C(:,:,Iz),faceAlpha);
    
end

function h=mybox(ax,XM,YM,ZM,C,facealpha)
% MYBOX: Draws surface on all 6 faces of the box defined by arguments
%
% USAGE:
%   h=mybox(XM,YM,ZM,C);
%
% TO 110501


h=NaN(6,1);

if length(size(squeeze(XM)))<3, return; end  % all matrices must be 3D to allow drawing

h(1)=mysurf(ax,XM(  1,:,:),YM(  1,:,:),ZM(  1,:,:),C(  1,:,:),facealpha);
h(2)=mysurf(ax,XM(end,:,:),YM(end,:,:),ZM(end,:,:),C(end,:,:),facealpha);
h(3)=mysurf(ax,XM(:,  1,:),YM(:,  1,:),ZM(:,  1,:),C(:,  1,:),facealpha);
h(4)=mysurf(ax,XM(:,end,:),YM(:,end,:),ZM(:,end,:),C(:,end,:),facealpha);
h(5)=mysurf(ax,XM(:,:,  1),YM(:,:,  1),ZM(:,:,  1),C(:,:,  1),facealpha);
h(6)=mysurf(ax,XM(:,:,end),YM(:,:,end),ZM(:,:,end),C(:,:,end),facealpha);

end

function h=mysurf(ax,XM,YM,ZM,C,facealpha)
% MYSURF: Plots a single surface defined by arguments
%
% USAGE:
%   h=mysurf(XM,YM,ZM,C)
%
% TO 110501

h=surf(ax,...
    squeeze(XM),...
    squeeze(YM),...
    squeeze(ZM),...
    squeeze(C ),...
    'facecolor','interp',...
    'edgecolor','none','facealpha',facealpha);

end
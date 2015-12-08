function h=plotMeshUpdate(o,h,Var,varargin)
% h = gr.plotMeshUpdate(h,Var,varargin)
% Update the mesh of plotMesh with new data Var
%
% TO 120619

if nargin<2, error('%s: handles required as second input',mfilename); end
if nargin<3, Var = o.ZM; end

if ~all(size(Var) == o.size),
    error('%s: Size of data [%d %d %d] must comply with size of grid [%d %d %d]',mfilename,size(Var),gr.size);
end

if nargin<4, varargin={'visible','on'}; end
    
if exist('var','var')
    ix=1;      set(h(1),'Cdata',squeeze(Var(:,ix,:)),varargin{:});
    ix=o.Nx+1; set(h(2),'Cdata',squeeze(Var(:,max(1,ix-1),:)),varargin{:});
    iy=o.Ny+1; set(h(3),'Cdata',squeeze(Var(max(1,iy-1),:,:)),varargin{:});
    iy=1;      set(h(4),'Cdata',squeeze(Var(iy,:,:)),varargin{:});
    iz=o.Nz+1; set(h(5),'Cdata',squeeze(Var(:,:,max(1,iz-1))),varargin{:});
    iz=1;      set(h(6),'Cdata',squeeze(Var(:,:,iz)),varargin{:});
end

drawnow;

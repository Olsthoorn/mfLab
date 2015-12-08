function plotXS(gr,varargin)
% gr.plotXS(ax,lineSpec,[iy[,varargin]]); % plot the grid of vertical section through iy.
% 
% TO 130306

if nargin<2
    lineSpec = 'c';
elseif ~isLineSpec(varargin{1})
    lineSpec = 'c';
else
    lineSpec = varargin{1};
    varargin(1)=[];
end

if ~isempty(varargin) && isaxis(varargin{1})
    ax = varargin{1}; varargin(1)=[];
else
    ax = gca;
end

set(ax,'nextplot','add');

if isempty(varargin)
    iy=1;
else
    if isnumeric(varargin{1})
        iy=varargin{1}(1);
        varargin(1)=[];
    else
        iy=1;
    end
end
   
if iy<1 || iy>gr.Ny
    error('%s: iy=%d (rowNr to plot XSec of must be between 1 and Ny=%d',...
        mfilename,iy,gr.Ny);
end

xGr = [gr.xGr(1:end-1); gr.xGr(2:end)]'; xGr=xGr(:)';
zxs = XS(gr.Z(iy,:,:));

for iz = 1:size(gr.Z,3)
    z =zxs([iz iz],:)'; z=z(:)';
    plot(ax,xGr,z,lineSpec,varargin{:});
end

for ix=1:gr.Nx
    plot(ax,xGr([ix+0 ix+0]),zxs([1 end],ix),lineSpec,varargin{:});
    plot(ax,xGr([ix+1 ix+1]),zxs([1 end],ix),lineSpec,varargin{:});
end



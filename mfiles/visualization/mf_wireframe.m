function mf_wireframe(ax,xGr,yGr,zGr,clr)
%MF_WIREFRAME plots wireframe definend by arguments
%
% Exmple:
%   mf_wireframe(xGr,yGr,zGr[,color])
%
%  one value for each coordinate must be used to determine
%  position in space
%
% ToDo: give proper example
%
% See also: mf_3Dblock
%
% TO 110501

set(ax,'nextplot','add');

grey=[0.8 0.8 0.8]; % use true color

if nargin<4, clr=grey; end

if nargin<3 || ~exist('zGr','var'), zGr=0; end


%% Analyze all possible forms of input

if any(size(xGr)==numel(xGr)) % then xGr is vector
    Nx=numel(xGr);
else
    Nx=size(xGr,2); % then xGr must be 2D or 3D, assume right shape
end

if any(size(yGr)==numel(yGr)),  % yGr is vector
    Ny=numel(yGr);
else
    Ny=size(yGr,1); % then yGr must be 2D or 3D array
    if size(yGr,2) ~=Nx
        error('mfLab:mf_wireframe:inputdim','input dimension yGr does not meet that of xGr');
    end
end

if any(size(zGr)==numel(zGr)) % then zGr is a vector
    Nz=numel(zGr);
else
    Nz=numel(zGr)/(Ny*Nx); % else, zGr ia a 2D or 3D array
    if size(zGr,1)~=Ny || size(zGr,2) ~=Nx,
        error('mfLab:mf_wireframe:inputdim','input dimension zGr does not match xGr and or yGr');
    end
end

%% Now generate full size 3D arrays for all coordinates

if numel(xGr)<Ny*Nx
    xGr=repmat(xGr(:)',[Ny,1,Nz]);
elseif numel(xGr)==Nx*Ny
    xGr=repmat(xGr,[1,1,Nz]);
end

if numel(yGr)<Nx*Ny
    yGr=repmat(yGr(:),[1,Nx,Nz]);
elseif numel(yGr)==Nx*Ny
    yGr=repmat(yGr,[1,1,Nz]);
end

if numel(zGr)<Nx*Ny*Nz
    if any(size(zGr)==numel(zGr))
        zGr=reshape(zGr(:),[1,1,numel(zGr)]);
        zGr=repmat(zGr,[Ny,Nx,1]);
    end
end    

%% Finally plot the grid, so that the coordinates can be anything, regular
% or irregular

for ix=1:size(xGr,2)
    for iy=1:size(yGr,1)
        for iz=1:size(zGr,3)-1
            line(...
                squeeze(xGr(iy,ix,iz+[0 1])),...
                squeeze(yGr(iy,ix,iz+[0 1])),...
                squeeze(zGr(iy,ix,iz+[0 1])),...
                'color',clr,'parent',ax);
        end
    end
end

for iy=1:size(yGr,1)
    for iz=1:size(zGr,3)
        for ix=1:size(xGr,2)-1
            line(xGr(iy,ix+[0 1],iz),yGr(iy,ix+[0 1],iz),zGr(iy,ix+[0 1],iz),'color',clr,'parent',ax);
        end
    end
end

for iz=1:size(zGr,3)
    for ix=1:size(xGr,2)
        for iy=1:size(yGr,1)-1
            line(xGr(iy+[0 1],ix,iz),yGr(iy+[0 1],ix,iz),zGr(iy+[0 1],ix,iz),'color',clr,'parent',ax);
        end
    end
end

drawnow;
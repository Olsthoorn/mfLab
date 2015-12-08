function plotGrid(o,varargin)
%% h = o.plotGrid(ax,clr,well,z,figname,figcoords)
% PLOTGRID: Plots the grid lines in color clr given the coordinates xGr yGr
%    and possibly well locations in blue if well is a struct or wellObj
%    whose elements have fields x and y
%
% 'zDir',true
% USAGE:
%   griDObj.plotGrid(plotOptions) % plot xy grid in 2D with cyan lines
%   gridObj.plotGrid(ax,..
%   gridObj.plotGrid(..., 'x',x,...    plots zy grid at x=x
%   gridObj.plotGrid(..., 'y',y,...    plots zx grid at y=y
%   gridObj.plotGrid(..., 'z',z,...    plots xy grid at z=z
%   gridObj.plotGrid(..., 'ix',ix,...  plots zy grid at ix=ix
%   gridObj.plotGrid(..., 'iy',iy,...  plots zx grid at iy=iy
%   gridObj.plotGrid(..., 'iz',iz,...  plots xy grid at iz=iz
%   gridObj.plotGrid(..., wellObj, ... also plots wells
%   gridObj.plotGrid(..., property,value ...
%   gridObj.plotGrid(..., 'figPos' ,figPos, ...
%   gridObj.plotGrid(..., 'xLim', xLim, ... 'yLim',yLim, ... 'zlim',zlim ...
%    additonally add graphic properties at the end of the input like:
%   gridObj.plotGrid(...,  'color','k','edgealpha','0.2','lineStyle','--');
%
% the x,y,z,ix,iy and iz are obligatory, they are used to select the layer
% row or column to be plotted as coordinates or index and also to select
% which wells to plot in case wellObj are specified. These wells are only
% plotted when the fall within the column row or layer for which the grid
% is plotted.
%
% SEE ALSO: well.plot3D gridObj.googleMap plotgrid
% TO 091201 100115 121016 151125

[xlim, varargin] = getProp(varargin,'xlim',[]);
[ylim, varargin] = getProp(varargin,'ylim',[]);
[zlim, varargin] = getProp(varargin,'zlim',[]);
[xLoc, varargin] = getProp(varargin,'x',[]);
[yLoc, varargin] = getProp(varargin,'y',[]);
[zLoc, varargin] = getProp(varargin,'z',[]);
[  ix, varargin] = getProp(varargin,'ix',[]);
[  iy, varargin] = getProp(varargin,'iy',[]);
[  iz, varargin] = getProp(varargin,'iz',[]);
[well, varargin] = getType(varargin,'wellObj',[]);
[well, varargin] = getType(varargin,'MNW1Obj',well);
[well, varargin] = getType(varargin,'MNW2Obj',well);
[ax, varargin]   = getType(varargin,'axis',[]);
[ax, varargin]   = getProp(varargin,'axis',ax);
[pos,varargin]   = getProp(varargin,'figPos',[]);

axProps = {'nextplot','add','xgrid','on','ygrid','on','zgrid','on'};

if isempty(ax)
    if isempty(pos)
        figure;
    else
        figure('position',pos);
    end
    ax = axes(axProps{:});
end

if rem(numel(varargin),2)==1
    [lineSpec, varargin] = getNext(varargin,'char','c');
    [IL,c,~,LS] = isLineSpec(lineSpec);
    if IL
        lineSpec = [c LS];
    else        
        lineSpec = 'b-';
    end
else
    lineSpec = 'b-';
end

%% plot zy

if ~isempty(xLoc) % plot zy section
    ix = hit(o.xGr,xLoc);
    plotzy(ax,ix);
    
    if ~isempty(well)
        for iw = numel(well):-1:1
            I(iw) = any(well(iw).Ix == ix);
        end
        if ~isempty(I)
            well(I).plotYS(ax);
        end
    end
    
    return
end

%% plot zx

if ~isempty(yLoc) % plot zx section
    iy = hit(o.yGr,yLoc);
    plotzx(ax,iy);
    
    if ~isempty(well)
        for iw = numel(well):-1:1
            I(iw) = any(well(iw).Iy == iy);
        end
        well(I).plotXS(ax);
    end
    
    return;
end

%% plot yx

if ~isempty(zLoc)
    plotyx(ax,zLoc);

    if ~isempty(well)
        for iw=numel(well):-1:1
            I(iw) = well(iw).iz(1)>=zLoc & well(iw).iz(end)<zLoc;
        end
        well(I).plotXY(ax);
    end
    
    return;
end
    
%% plot yz

if ~isempty(ix) % plot zy section
    plotzy(ax,ix);
    
    if ~isempty(well)
        for iw = numel(well):-1:1
            I(iw) = any(well(iw).Ix == ix);
        end
        well(I).plotYS(ax);
    end
    
    return;
end

%% plot xz

if ~isempty(iy) % plot zx section
    plotzx(ax,iy);
    if ~isempty(well)
        for iw = numel(well):-1:1
            I(iw) = any(well(iw).Iy == iy);
        end
        well(I).plotXS(ax);
    end
    return;
end

%% plot yx

if ~isempty(iz) % plot xy section
    plotyx(ax);
    if ~isempty(well)
        for iw = numel(well):-1:1
            I(iw) = any(well(iw).Iz == iz);
        end
        well(I).plotXY(ax);
    end
    return
end
   
%% functions

function plotyx(ax)
    if isempty(xlim)
        xlim = o.xGr([1 end]);
    end
    if isempty(ylim)        
        ylim = o.yGr([end 1]);
    end
        set(ax,'xlim',xlim,'ylim',ylim);
    for ix_ = 1:o.Nx+1
        plot(ax,o.xGr([ix_ ix_]),o.yGr([1 end]),lineSpec,varargin{:});
    end
    for iy_ = 1:o.Ny+1
        plot(ax,o.yGr([iy_ iy_]),o.xGr([1 end]),lineSpec,varargin{:});
    end
end
    
function plotzy(ax,ix)
    if isempty(xlim)
        xlim = o.yGr([end 1]);
    end
    if isempty(zlim)
        zlim = o.yGr(([end 1]));
    end
    set(ax,'xlim',xlim,'ylim',zlim);
    zGr = zeros(2,o.Ny);
    yGr = [o.yGr; o.yGr];
    for iz_ = 1:o.Nz+1
        zGr(1,:) = o.Z(:,ix,iz_);
        zGr(1,:) = o.Z(:,ix,iz_);
        plot(ax,yGr(2:end-1),zGr(:)',lineSpec,varargin{:});
    end
    zTop = o.Z(:,ix,  1); zGrt = max([zTop([1 1:end]); zTop([1:end end])]);
    zBot = o.Z(:,ix,end); zGrb = min([zBot([1 1:end]); zBot([1:end end])]);
    for iy_ = 1:o.Ny+1
        plot(ax,o.yGr([iy_ iy_]),[zGrb(iy_) zGrt(iy_)],lineSpec,varargin{:});
    end
end
function plotzx(ax,iy)
    if isempty(xlim)
        xlim = o.xGr([1 end]);
        set(ax,'xlim',xlim);
    end
    if isempty(zlim)
        zlim = ([min(o.zGr(:)) max(o.zGr(:))]);
    end
    set(gca,'xlim',xlim,'ylim',zlim);
    
    zGr = zeros(2,o.Nx);
    xGr = [o.xGr; o.xGr];
    for iz_ = 1:o.Nz+1
        zGr(1,:) = o.Z(iy,:,iz_);
        zGr(2,:) = o.Z(iy,:,iz_);
        plot(ax,xGr(2:end-1),zGr(:)',lineSpec,varargin{:});
    end
    zTop = o.Z(iy,:,  1); zGrt = max([zTop([1 1:end]); zTop([1:end end])]);
    zBot = o.Z(iy,:,end); zGrb = min([zBot([1 1:end]); zBot([1:end end])]);
    for ix_ = 1:o.Nx+1
        plot(ax,o.xGr([ix_ ix_]),[zGrb(ix_) zGrt(ix_)],lineSpec,varargin{:});
    end
end

end

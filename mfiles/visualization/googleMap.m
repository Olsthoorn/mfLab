function h = googleMap(varargin)
%GOOGLEMAP  geneates a google Map object
%
% Example:
%    GM = gr.googleMap(ax,xLLUR ,yLLUR,z,mapType [,coordSystem]);
%    GM = gr.googleMap(E_LLUR,N_LLUR,z,mapType [,coordSystem]);
%
% INPUT
%    x = coordinates of LL and UR of picture;
%    y = coordinates of LL and UR of picture;
%    z = elevation of picture [m];
%    maptype     = [{roadmap} satellite hybrid terrain]
%    coordSystem = [ {wgs} rd]
%        rd = Dutch National Rijksdriehoekssystem 
%
%     Notice that coordinate must be N and E (lat long) when wgs is used
%     and in m when rd is used.
%
% ToDo: make proper example
%
% TO 121016 130402

MAXZOOM = 21;
MAXPIX  = 640;

[ax,varargin] = getNext(varargin,'axis',[]);
[ax,varargin] = getProp(varargin,'axis',ax);
if isempty(ax)
    ax = gca;
end

[x, varargin] = getNext(varargin,'double',[]);

if isempty(x) || numel(x) ~=2 || abs(diff(x))<50*eps
    error('%s: first argument x must be two distinct x-coordinates',mfilename);
end

[y, varargin] = getNext(varargin,'double',[]);

if isempty(y) || numel(y) ~=2 || abs(diff(y))<50*eps
    error('%s: second argument y must be two distinct y-coordinates',mfilename);
end

[z, varargin] = getNext(varargin,'double',[]);

if isempty(z) || numel(z) ~=1
    error('%s: thrid argument z must be a single z value (elevation of picture)',mfilename);
end

[mapType    ,varargin] = getNext(varargin,'char','roadmap');
[coordSystem,varargin] = getNext(varargin,'char','wgs');
% exchange if the user choose the wrong order
if strmatchi(coordSystem,{'roadmap','satellite','terrain','hybrid'})
    dum         = coordSystem;
    coordSystem = mapType;
    mapType     = dum;
end
    

switch lower(coordSystem)
    case 'rd'
        [rdLL(2) rdLL(1)] = rd2wgs(min(x),min(y));
        [rdUR(2) rdUR(1)] = rd2wgs(max(x),max(y));
    case {'wgs','wgs84'}
        rdLL = [min(y),min(x)];
        rdUR = [max(y),max(x)];
        
    otherwise
        error('%s: Unknown coordinate system <<%s>>',mfilename,coordSystem);
end
        
center  = 0.5*(rdLL+rdUR);

gpLL = googlePointObj(rdLL(1),rdLL(2),MAXZOOM);
gpUR = googlePointObj(rdUR(1),rdUR(2),MAXZOOM);

for iZoom = MAXZOOM:-1:0
    gpLL = googlePointObj(gpLL.Lat,gpLL.Lon,iZoom);
    gpUR = googlePointObj(gpUR.Lat,gpUR.Lon,iZoom);
    Nx = (gpUR.ix*255+gpUR.px) - (gpLL.ix*255+gpLL.px);
    Ny = (gpLL.iy*255+gpLL.py) - (gpUR.iy*255+gpUR.py);
    if Nx<=MAXPIX && Ny<=MAXPIX,
        break;
    end
end

GM = googleMapObj(center,iZoom,[Nx Ny],mapType,'PNG');

xp = linspace(min(x),max(x),Nx+1);
yp = linspace(max(y),min(y),Ny+1);

[rgb, ~       ] =  getNext(varargin,'char','rgb');

hold(ax,'on');

if strcmpi(rgb,'rgb')
    set(get(ax,'parent'),'renderer','zbuffer');
    h = surf(ax,xp,yp,z*ones(Ny+1,Nx+1),map2rgb(GM.X,GM.iminfo.Colormap),'edgecolor','none','parent',ax);
else
    set(get(ax,'parent'),'renderer','painters'); % guarantees flicker free video
    h = surf(ax,xp,yp,z*ones(Ny+1,Nx+1),double(GM.X),'edgecolor','none','parent',ax);
    colormap(GM.iminfo.Colormap);
end
shading('flat');

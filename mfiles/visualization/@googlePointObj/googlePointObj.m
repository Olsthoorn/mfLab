classdef googlePointObj
% DESCRIPTION:
%   googlePointObj holds the coordinates of points
%   that consist of 5 number:
%      zoom level, tile coordinates (ix,iy) and pixel coordinates (px,py) withing
%      the tile.
%    See:
%      http://maps.google.com/maps/api/staticmap
%
% SEE ALSO: googleMapObj
%
% TO 110501 121005 121008

    properties
        Lat,Lon
        zoom
        name
        Nx = 256;
        Ny = 256;
    end
    properties (Dependent = true)
        w
        phi,lam
        x,y
        ix,iy
        px,py
    end
    methods
        function o=googlePointObj(Lat,Lon,zoom)
            % GP = googlePointObj(Lat,Lon,zoomLevel)
            % get a google Map point object
            % TO 121016
            
            if nargin==0
                return;
            end
            if numel(Lat) ~= numel(Lon)
                error('%s: numel of Lat (%d) ~= numel of Lon (%d)', ...
                    mfilename,numel(lat),numel(Lon));
            end
            if numel(Lat) ~= numel(zoom) && numel(zoom) ~=1
                error('%s: numel(Lat) and numel(Lon)=%d must equal numel(zoom)=%d or numel(zoom) must be 1',...
                    mfilename,numel(Lon),numel(zoom));
            end
            
            for i=length(Lat):-1:1
                o(i).Lat = Lat(i);
                o(i).Lon = Lon(i);
                o(i).zoom= zoom(min(numel(zoom),i));
            end
        end
        
        function w = get.w(o),     w =2^-o.zoom; end
        function phi = get.phi(o), phi = pi*o.Lat/180; end % in vert plane
        function lam = get.lam(o), lam = pi/180*o.Lon; end % angle hor plane
        function y   = get.y(o)
            y   = 0.5*(1-log((1+sin(o.phi))./(1-sin(o.phi)))/(2*pi)); end
        function x   = get.x(o),   x   = (o.lam/pi+1)/2; end
        function ix  = get.ix(o),  ix  = fix(o.x/o.w); end
        function iy  = get.iy(o),  iy  = fix(o.y/o.w); end
        function px  = get.px(o),  px =round(255*(o.x/o.w-o.ix)); end
        function py  = get.py(o),  py =round(255*(o.y/o.w-o.iy)); end
    end
end

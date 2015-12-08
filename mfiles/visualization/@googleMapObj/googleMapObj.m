classdef googleMapObj
% DESCRIPTION:
%   Google Map images may be retrieved through the Google Maps API for
%   static maps. This is done by constructing the correct URL and dropping
%   it into the URL bar of your browser.
%
%    Google maps image requests have the following form:
%    See:
%      http://maps.google.com/maps/api/staticmap?parameters
%
% Google API site example
% URL=[...
%     'http://maps.google.com/maps/api/staticmap?',...
%     'center','Brooklyn+Bridge,New+York,NY','&',...
%     'zoom=','14','&',...
%     'size=','512x512','&',...
%     'maptype=','roadmap','&',...
%     'markers=','color:blue%7Clabel:S%7C40.702147,-74.015794','&',...
%     'markers=','color:green%7Clabel:G%7C40.711614,-74.012318','&',...
%     'markers=','color:red%7Ccolor:red%7Clabel:C%7C40.718217,-73.998284','&',...
%     'sensor=''false'];
%
% CORRESPONDING MATLAB STRUCT url
%    see selftest code how it is made in a Matlab struct
%
% SE ALSO: mf_GM2PNG, mf_GMTiles mf_GMLL2pix kmlpath mf_kmlpath2RD wgs2rd rd2wgs  
%
% TO 110501 121005

    properties
        R2      = 6371007.2;  % The Earth's authentic radius (Geodetic Union);
        name    = '',
        center  = 'Amsterdam',  % Lat Lon or name of place
        zoom    = 16,
        pixels  = [640 640], % picture size in pixels
        maptype = 'satellite',
        format  = 'png',
        path    = [],  % struct
        paths   = '',  % URL part holding eventual paths
        marker  = [],  % struct
        markers = '',  % URL part holding eventual markers
        URL     = '',  % URL part holding URL
        X     % the image als 2D matrix
        A     % the image as RGB
        iminfo  % iminfo about the image
        map
    end
    properties (Dependent=true)
        w         % relative width and height of tile 2^(-o.zoom)
        gpc       % googleMapPoint of center
        gpLL,gpUR % googleMapPoints of LL and UR of tile
        tile      % tile properties witdth and height dx dy
    end
    methods
        function o=googleMapObj(center,zoom,pixels,maptype,format,marker,path)
            if nargin<3
                 return;
            end
            
            if ischar(center) % the argument is the name of a kml file
                try
                    o.name = center;
                    [Lat Lon]=kmlpath(o.name); % LON LAT
                    o.center = [Lat Lon];     % LAT LON
                catch ME %#ok
                    o.name= center;
                    [Lat Lon] = geoAddress(o.name);
                    o.center = [Lat Lon];
                end
            elseif isstruct(center)
                o.name = center;
                [Lat Lon o.name] = geoAdress(o.name);
                o.center = [Lat Lon];
            elseif size(center,1)==2
                o.center = mean(center,1);
            elseif isnumeric(center)
                o.center = center;  % LAT LON
            else
                error('%s: center has illegal format',fmilename);
            end
            % center may also be a kml file representing a single pin
            o.zoom   = min(max(0,zoom),21);
            o.pixels = [min(640,pixels(1)), min(640,pixels(end))];
            %% Verify maptype input
            if nargin>=4
                o.maptype=maptype;
                if ~ismember(o.maptype,{'roadmap','satellite','terrain','hybrid'})
                    error('Illegal maptype %s\n',o.maptype);
                end
            end
            
            %% Picture format
            if nargin>=5
                o.format=lower(format);
                if ~ismember(o.format,{'png8','png','png32','gif','jpg','jpg-baseline'})
                    error('Illegal file format %s\n',o.format);
                end
            end
            
            if nargin>=6,
                o=setMarkers(o,marker);
            end
            
            if nargin>=7
                o=setPaths(o,path);
            end
                 
            % Necessary requests to obtain a sufficiently detailed picture            
            % assemble image file from pieces
            % Under construction
            o=o.setURL();        
            o=o.getImage();
        end
                
       function o=setURL(o)
           
            o.URL=[...
                'http://maps.googleapis.com/maps/api/staticmap?',...
                'center=',sprintf('%.6f,%.6f',o.center),'&',...
                'zoom=',sprintf('%d',o.zoom),'&',...
                'size=',sprintf('%dx%d',o.pixels),'&',...
                'format=',o.format,'&',...
                'maptype=',o.maptype,'&',...
                 o.paths,...
                 o.markers,...
                'sensor=false'];
            o.URL(o.URL==' ')='+';
       end
        
        function w    = get.w(o), w=2^(-o.zoom); end
        function gpc  = get.gpc( o)
            gpc =  googlePointObj(o.center(1),o.center(2),o.zoom);
        end
        function gpLL = get.gpLL(o)
            gpLL = o.gpc.movePix(-o.pixels(1)/2,+o.pixels(2)/2);
        end
        function gpUR = get.gpUR(o)
            gpUR = o.gpc.movePix( o.pixels(1)/2,-o.pixels(2)/2);
        end            
        function o=web(o)
            % mapObj.web -- show image on Matlab's browser
            web(o.URL);
        end
        function tile = get.tile(o)
            % size of the figure in meters (not the google tile)
            tile.Lx = o.R2 * (o.gpUR.lam-o.gpLL.lam) * cos(o.gpc.phi);
            tile.Ly = o.R2 * (o.gpUR.phi-o.gpLL.phi);
            tile.dx = tile.Lx/256;
            tile.dy = tile.Ly/256;
        end

        function o=browser(o)
            % mapObj.browser() -- show image on system browser
            web(o.URL,'-browser');
        end
        function o=setMarkers(o,marker)
        % mapOb = mapObj.setMarkers(markers);
        % URLMARKERS: Generate URL string part with markers
        % also store marker struct (with fields color,label,LATLON)
        %
        % USAGE:
        %   o=o.setMarkers(url)
        %
        % TO 110501 111004
        
            if nargin>1
                o.marker=marker; % o.marker = strucct, o.markers is urlstring
            end     
            
            if ~isempty(o.marker)
                o.markers='';
                for i=1:length(marker)
                    o.markers=[o.markers,...
                     'markers=',...
                     'color:',marker(i).color,'%7C',...
                     'label:',marker(i).label,...
                      sprintf('%%7C%.6f,%.6f',marker(i).LATLON),...
                      '&'];  % included so that if markers is not a field no crashes occur
                end
            end
        end
        
       function o=addPath(o,kmlfilename,linecolor,weight,fillcolor)
            if isempty(o.path), i=1; else i=length(o.path)+1; end
            o.path(i).fname=kmlfilename;
            o.path(i).linecolor=linecolor;
            o.path(i).weight = weight;
            o.path(i).fillcolor = fillcolor;
            o.path(i).LATLON     =  kmlpath(o.path(i).fname);
            
            o.paths=o.setPath(o.path);
        end

        function o=setPaths(o,path)
            % adds string to URL that define a path that can be plotted
            %
            % path must be a struct array with the following fields
            %   path.color
            %   path.weightstyles
            %   path.locations = [lat long; lat long; et]
            % TO 110501
        
            if nargin>1
                o.path=path;
            end
            
            if ~isempty(o.path)
                o.paths='';
                for i=1:length(o.path)
                    if isfield(o.path(i),'color')
                        color=sprintf('color:%s',o.path(i).color);
                    else
                        color='';
                    end
                    if ~isempty(o.path(i).weight)
                        weight=sprintf('%%7Cweight:%g',o.path(i).weight);
                    else
                        weight='';
                    end
                    if ~isempty(o.path(i).fillcolor)
                        fillcolor=sprintf('%%7Cfillcolor:%s',o.path(i).fillcolor);
                    else
                        fillcolor='';
                    end
                    o.paths=[o.paths,'path=',...
                        color,...
                        weight,...
                        fillcolor,...
                        sprintf('%%7C%.6f,%.6f',path(i).LATLON'),...
                        '&'...
                      ];
                end          
            end
        end

        function o = image(o,xlim,ylim)
            % GM.image(xlim,ylim)
            % GM.image();
            % show the image
            % TO 121006
            figure; hold on;
            colormap(o.iminfo.Colormap);
            if nargin<2
                image(o.X);
                xlabel('px [pixels]'); ylabel('py [pixels]');
            else
                image(xlim,ylim,o.X(end:-1:1,:,:));
                axis('tight');
                set(gca,'ydir','normal');
                xlabel('x [m]'); ylabel('y[m]');
            end
            title(sprintf('%s, zoomlevel = %d',o.name,o.zoom));
        end

        function o = getImage(o)
            o.A      = imread( o.URL,o.format);
            [o.X,o.map]= imread( o.URL,o.format);
            o.iminfo = imfinfo(o.URL,o.format);
        end
        
        function o = getMultiImage(o)
            % problem: how to get multiple figures without a copyright?
            [o.X,o.map] = imread(o.URL,oformat);
            o.iminfo = imfinfo(o.URL,o.format);
        end
        
        function o=plotRDImage(o)
            [xLL,yLL]=wgs2rd(o.LL(2),o.LL(1));
            [xUR,yUR]=wgs2rd(o.UR(2),o.UR(1));

            image([xLL xUR],[yLL yUR],o.A);
            set(gca,'ydir','normal');
           
        end
    end
end

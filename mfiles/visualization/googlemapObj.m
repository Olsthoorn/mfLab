classdef googlemapObj
%GOOGLEMAPOBJ class def for google map objects
%
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
% TO 110501

    properties
        name    = '',
        center  = 'Amsterdam',  % Lat Lon or name of place
        zoom    = 16,
        size    = [640 640],
        maptype = 'satellite',
        format  = 'png',
        path    = [],  % struct
        paths   = '',  % URL part holding eventual paths
        marker  = [],  % struct
        markers = '',  % URL part holding eventual markers
        URL     = '',  % URL part holding URL
        xLim, yLim, dx, dy,
        LL = NaN(1,2); % Lower left of image in Lat Lon
        UR = NaN(1,2); % Upper right of image in Lat Lon
        X     % the image als 2D matrix
        A     % the image as RGB
        iminfo  % iminfo about the image
    end
    methods
        function o=googlemapObj(center,zoom,size,maptype,format,marker,path)
            if nargin<3
                o=o.selftest();
                return;
            end
                        
            %% Correct form of center [LAT LON] or character location
            % Check zoom and size input
            o.center = center;  % LAT LON
            o.zoom   = min(max(0,zoom),21);
            o.size   = [min(640,size(1)), min(640,size(end))];
            
            %% Verify maptype input
            if nargin>=4
                o.maptype=maptype;
                if ~ismember(o.maptype,{'roadmap','satellite','terrain','hybrid'})
                    error('Illegal maptype %s\n',o.maptype);
                end
            end
            
            %% Picture format
            if nargin>=5
                o.format=format;
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
            
            o=o.setBox();
            
            o=o.setURL();
            
            o=o.getImage();
        end
                
       function o=setURL(o)
           
            o.URL=[...
                'http://maps.googleapis.com/maps/api/staticmap?',...
                'center=',sprintf('%.6f,%.6f',o.center),'&',...
                'zoom=',sprintf('%d',o.zoom),'&',...
                'size=',sprintf('%dx%d',o.size),'&',...
                'format=',o.format,'&',...
                'maptype=',o.maptype,'&',...
                 o.paths,...
                 o.markers,...
                'sensor=false'];
            o.URL(o.URL==' ')='+';
       end
        
        function o=web(o)
            % mapObj.web -- show image on Matlab's browser
            web(o.URL);
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
                    o.markes=[o.markers,...
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
        % URLPATH: Generate URL string part with path
        %
        % USAGE:
        %   o=o.setpath(path)
        %
        % path must have
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
                    o.paths=[o.path,'path=',...
                        color,...
                        weight,...
                        fillcolor,...
                        sprintf('%%7C%.6f,%.6f',path(i).LATLON'),...
                        '&'...
                      ];
                end          
            end
        end
        
        function o=setBox(o)
            % setUR and LL of object in LL from center, zoom and size
            
            [px,py,ix,iy]=mf_GMLL2pix(o.center(1),o.center(2),o.zoom);

            ixLL=ix; iyLL=iy; ixUR=ix; iyUR=iy;

            pxLL=px-o.size(1)/2; while pxLL<  0, ixLL=ixLL-1; pxLL=pxLL+256; end
            pxUR=px+o.size(1)/2; while pxUR>256, ixUR=ixUR+1; pxUR=pxUR-256; end 
            pyLL=py+o.size(2)/2; while pyLL>256, iyLL=iyLL+1; pyLL=pyLL-256; end
            pyUR=py-o.size(2)/2; while pyUR<  0, iyUR=iyUR-1; pyUR=pyUR+256; end

            [o.LL(1) o.LL(2)]=mf_GMpix2LL(ixLL,iyLL,o.zoom,pxLL,pyLL);
            [o.UR(1) o.UR(2)]=mf_GMpix2LL(ixUR,iyUR,o.zoom,pxUR,pyUR);
        end
        
        function o=getImage(o)
            o.X=imread(o.URL,o.format);
            o.X=o.X(end:-1:1,:);
            o.iminfo=imfinfo(o.URL);
            o.A=ind2rgb(o.X,o.iminfo.Colormap);
        end
        
        function o=plotRDImage(o)
            [xLL,yLL]=wgs2rd(o.LL(2),o.LL(1));
            [xUR,yUR]=wgs2rd(o.UR(2),o.UR(1));

            image([xLL xUR],[yLL yUR],o.A);
            set(gca,'ydir','normal');
           
        end

        function o=selftest(o)
            % SELFTEST: Runs self test if mf_GM2PNG is called without arguments
            %
            % USAGE:
            %   selftest
            %
            % TO 110501

            mrkr(3).color='red';
            mrkr(3).label='C';
            mrkr(3).LATLON=[40.718217,-73.998284'];

            mrkr(2).color='green';
            mrkr(2).label='G';
            mrkr(2).LATLON=[40.711614,-74.012318];

            mrkr(1).center='Brooklyn+Bridge,New+York,NY';
            mrkr(1).zoom=14;
            mrkr(1).size=[512 512];
            mrkr(1).maptype='roadmap';

            mrkr(1).color='blue';
            mrkr(1).label='S';
            mrkr(1).LATLON=[40.702147,-74.015794];
            
            o=o.setMarkers(mrkr);

            o=o.setURL(o); % make url and show map

            fprintf('YOUR URL is\n%s\n',o.URL); % and show it
        end
        
        function o=mf_GM_tile_size(o)

            o.xLim=o.dx*(o.size(1)/256)*[-0.5 0.5];  % every GM tile is 256x256 pixels
            o.yLim=o.dy*(o.size(2)/256)*[-0.5 0.5];  % of known absolute size

            %% Here we use the URL twice. First to get the figure into Matlab
            % Then to get its information, becuase we need the colormap later on
            o.A      = imread(o.URL);                 % get it directly into matlab
            o.iminfo = imfinfo(o.URL);               % get associated image iminfo

            %% Here the fiture is plotted, using the coordianates in m centered around
            %% the center of the figure.
            image(o.xLim,o.yLim([2 1]),o.A); colormap(o.iminfo.Colormap); % Plot it using the colormap from the iminfo

            hold on;
            set(gca,'ydir','normal');
            xlabel('x [m relative to tile center]');
            ylabel('y [m relative to tile center]');
            title('A small Dubai World Island, suitable for subsurface water revervoirs?');

            %% Now lets use the GE path by converting its Lat Lon coordinates into x,y
            %% coordinates relative to the point url.center.
            %island=   mf_GMLL2XY(o.path(1).LATLON,o.center);
            %reservoir=mf_GMLL2XY(o.path(2).LATLON,o.center);

            %% We may add a square of 600 by 600 m to allow verification with Google
            %% Earth's ruler

            %square   =mf_GMLL2XY(o.path(3).LATLON,o.center);

            %% Here we show how this can be used in your model
        end      
    end
end

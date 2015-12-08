function [xW yS xE yN]=mf_GM_PIC_XY(url,varargin)
%MF_GM_PIC_XY retrieves a goolge map image and puts it in Lat Lon axis
%
% Example: to be worked out
%   mf_GM_PIC_XY
%    selftest + help
%
%   URL=mf_GM2PNG(url)  % same as URL=mf_GM2PNG(url) defaul PNG fomat
%     where url is a struct with the correct fields, see example
%     URL = provided
%
%   [xW yS xE yN]=mf_GM_PIC_XY(url,maptype);
%    as before, overwrites default maptype
%    maptype is one of'roadmap','satellite','terrain','hybrid';
%
%   [xW yS xE yN]=mf_GM_PIC_XY(url,maptype,figtype)
%    as before, overwrites als default figtype
%    figtype is one of 'png8','png','png32','gif','jpg','jpg-baseline'
%
% ISSUE: How to comply with Goole's Copyright? > User's responsability
%
% MOTIVATION:
%   Static Google Maps may be retrieve through the Google Maps APIL for static maps.
%   This is usefule for showing simulatons results on a correct map
%   background
%
% DESCRIPTION:
%   Google Map images may be retrieved through the Google Maps API for
%   static maps. This is done by constructing the correct URL and dropping
%   it into the URL locaiton of your browser.
%   Requests have the following form:
%
%   http://maps.google.com/maps/api/staticmap?parameters
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
% CORRESPONDING MALTAB STRUCT url
%    see selftest code how it is made in a Matlab struct
%
% SE ALSO: mf_GM2PNG, mf_GMTiles mf_GMLL2pix kmlpath mf_kmlpath2RD wgs2rd rd2wgs  
%
% TO 110501

if nargin==0, selftest; return; end

if ~isfield(url,'maptype'), url.maptype='satellite'; end % see legalmaptypes below
if ~isfield(url,'format'),  url.format ='png';       end % see lageformat below

if nargin>1
    url.maptype=varargin{1};
    legalmaptypes={'roadmap','satellite','terrain','hybrid'};
    try
        strmatchi(url.maptype,legalmaptypes);
    catch ME
        ME.message=sprintf('Second input arg must be one of %s\n',sprintf(' ''%s''',legalmaptypes{:}));
        error(ME);
    end
end
 
if nargin>2
    url.format=varargin{2};
    legalpictypes={'png8','png','png32','gif','jpg','jpg-baseline'};
    try
        strmatchi(url.format,legalpictypes);
    catch ME
        ME.message=sprintf('Third input arg must be one of %s\n',sprintf(' ''%s''',legalpictypes{:}));
        error(ME);
    end
end

URL=mf_GM2PIC(url,url.maptype,url.format);

[px,py,ix,iy]=mf_GMLL2pix(url.center(1),url.center(2),url.zoom);

pxW=px-url.size(1)/2;
pxE=px+url.size(1)/2;

pyN=py-url.size(2)/2;
pyS=py+url.size(2)/2;

xW = mf_GMpix2XY(ix,iy,url.zoom,pxW,py);
xE = mf_GMpix2XY(ix,iy,url.zoom,pxE,py);

[dum yS]=mf_GMpix2XY(ix,iy,url.zoom,px,pyS);
[dum yN]=mf_GMpix2XY(ix,iy,url.zoom,px,pyN);

dx=xE-xW; xW=-dx/2;  xE=+dx/2;
dy=yN-yS; yS=-dy/2;  yN=+dy/2;

A    = imread(URL);                % get it directly into matlab
info = imfinfo(URL);               % get associated image info

image([xW xE],[yN yS],A);
colormap(info.Colormap);           % its colormap from info
set(gca,'ydir','normal');
xlabel('X [m]');
ylabel('Y [m]');
hold on;

grid on
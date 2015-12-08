function URL=GM2PIC(url,varargin)

% GM2PNG: Retrieves a goolge map image as png file to show simulation results on map
%
% USAGE: to be worked out
%   mf_GM2PNG
%    selftest + help
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

URL=[...
    'http://maps.google.com/maps/api/staticmap?',...
    'center=',gmcenter(url.center),'&',...
    'zoom=',sprintf('%d',min(max(0,url.zoom),21)),'&',...
    'size=',sprintf('%dx%d',min(640,url.size(1)),min(640,url.size(end))),'&',...
    'format=',url.format,'&',...
    'maptype=',url.maptype,'&',...
    urlpath(url),...
    urlmarkers(url),...
    'sensor=false'];

URL(URL==' ')='+';

if nargout==0
    web(URL); %,'-browser');
else
    fprintf('YOUR URL is:\n%s\n',URL);
end

end

function uc=gmcenter(uc)
% GMCENTER: Makes sure center is in correct string form no matter how specified
%
% USAGE: uc=gmcenter(uc)
%    uc has following form
%       'address string like: 'Brooklyn+Bridge,New+York,NY'
%     or LATLON value pair like [40.702147,-74.015794] in demcial degrees
%
% TO 110501

    if ischar(uc)
        uc(uc==' ')='+';
    else
       try
           uc=sprintf('%.6f,%.6f',uc);
       catch ME
            fprintf('mfLab:mf_GMPGN:gmcenter URL.center must be string or [LAT LON] numerical values in demical degrees');
            error(ME);
        end
    end
end

function s=urlmarkers(url)
% URLMARKERS: Generate URL string part with markers
%
% USAGE:
%   um=urlmarkers(url)
%
% TO 110501
 
    if isfield(url,'markers')
        s='';
        for i=1:length(url.markers)
            s=[s,...
             'markers=',...
             'color:',url.markers(i).color,'%7C',...
             'label:',url.markers(i).label,...
              sprintf('%%7C%.6f,%.6f',url.markers(i).LATLON),...
              '&'];  % included so that if markers is not a field no crashes occur
        end
    else
        s='';
    end

end

function s=urlpath(url)
% URLPATH: Generate URL string part with path
%
% USAGE:
%   s=urlpath(url)
%
% url.path must have
%   path.color
%   path.weightstyles
%   path.locations = [lat long; lat long; et]
% TO 110501
 
    if isfield(url,'path')
        
        s='';
        for i=1:length(url.path)

        
            if isfield(url.path(i),'color')
                color=sprintf('color:%s',url.path(i).color);
            else
                color='';
            end
            if isfield(url.path(i),'weight')
                weight=sprintf('%%7Cweight:%g',url.path(i).weight);
            else
                weight='';
            end
            if isfield(url.path(i),'fillcolor')
                fillcolor=sprintf('%%7Cfillcolor:%s',url.path(i).fillcolor);
            else
                fillcolor='';
            end


            s=[s,'path=',...
                color,...
                weight,...
                fillcolor,...
                sprintf('%%7C%.6f,%.6f',url.path(i).LATLON'),...
                '&'...
              ];
        end
        
    else
        s='';
    end

end

function selftest
% SELFTEST: Runs self test if mf_GM2PNG is called without arguments
%
% USAGE:
%   selftest
%
% TO 110501

    url.center='Brooklyn+Bridge,New+York,NY';
    url.zoom=14;
    url.size=[512 512];
    url.maptype='roadmap';
    
    url.markers(1).color='blue';
    url.markers(1).label='S';
    url.markers(1).LATLON=[40.702147,-74.015794];

    url.markers(2).color='green';
    url.markers(2).label='G';
    url.markers(2).LATLON=[40.711614,-74.012318];
    
    url.markers(3).color='red';
    url.markers(3).label='C';
    url.markers(3).LATLON=[40.718217,-73.998284'];

    mf_GM2PNG(url); % make url and show map
    
    eval('help mf_GM2PNG');           % show help info
    URL=mf_GM2PNG(url);               % generate the URL
    fprintf('YOUR URL is\n%s\n',URL); % and show it
  end

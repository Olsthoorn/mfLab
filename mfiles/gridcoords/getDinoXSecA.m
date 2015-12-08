function A=getDinoXSecA(varargin)
%GETDINOXSECA retrieve geo(hydro)logcal cross section from www.dinoloket.nl
%
% Examples:
%   mf_getDinoXSec(xRD,yDR,mdl)
%   mf_getDinoXSec(xRD,yDR)      % assumes model is 'DGM'
%   mf_getDinoXSec(xyRD,mdl)
%   mf_getDinoXSec(xyRD)         % assumes model is 'DGM'
%   mf_getDinoXSec(kmlpath,mdl)
%   mf_getDinoXSec(kmlpath)      % assumes model is 'DGM'
%   mf_getDinoXSec;           % this help
%
%   Possible models are
%   DGM REGIS GR FR DR OV GE NH ZH UT GE ZL NB LB
%
% Example1:
%   xRD=[106187 169360 249194]';
%   yRD=[488914 452815 468087]';
%   mdl='DGM';
%   mf_getDinoXSec(xRD,yDR,mdl)
%
% Example2:
%   xRD=[106187 169360 249194]';
%   yRD=[488914 452815 468087]';
%   xyRD=[xRD yRD]
%   mdl='REGIS';
%   mf_getDinoXSec(xDR,mdl);
%
% Example3:
%   mdl='REGIS';
%   mf_getDinoXSec(kmlpath,mdl);
%
% Hing
%   Directly use GE coordinates:
%   add a path in GE
%   export the path as .kml file, e.d. to mypath.kml in your current dirctory
%   then do the following
%   [xRD,yRD]=kmlpath2rd('mypath.kml')
%   getDinoXSec(xRd,yRd,'REGIS');
%
%   Once you have the immage in your browser, you may save the image as a
%   .png file in any place you want. This hold for the cross section and
%   its legend separately.
%
% See also: kmlpath wsg2rd rd2wsg
%
%   TO 110427


%% Defaults
switch nargin
    case 0
        xRD=[106187 169360 249194];
        yRD=[488914 452815 468087];
        getDinoXSecA(xRD,yRD,'DGM');
    case 1  % assume kmlpath given and model will be DGM
        kmlname=pclean(varargin{1});
        if exist(kmlname,'file')
            [xRD,yRD]=kmlpath2rd(kmlname);
        else 
            xRD=varargin{1}(:,1);
            yRD=varargin{1}(:,2);
        end
        mdl='DGM'; % first argument
        getDinoXSecA(xRD,yRD,mdl);
    case 2 % if ischar(varargin{1}) assume first arg is kmlpath name otherwise assume first arg is [xRd yRD]
        kmlname=pclean(varargin{1});
        if exist(kmlname,'file')
            [xRD,yRD]=kmlpath2rd(kmlname);
            mdl=varargin{2};
        else 
            xRD=varargin{1}(:,1);
            yRD=varargin{1}(:,2);
            mdl='DGM';
        end
    case 3
        xRD=varargin{1};
        yRD=varargin{2};
        mdl=varargin{3};
    otherwise
        eval(['help ' mfilename]);
        return
end

%% getting the xsection
switch mdl
    case 'DGM'
       modelName='LKN-DGM-v1.3';
       description='Landelijk model DGM v1.3 - 2009';
    case 'REGIS'
        modelName='lkn-regis-2.1';
        description='Landelijk model REGIS II.1 - 2008';
    otherwise
        modelName=['HGM-',mdl,'-HYP'];
        switch mdl
            case 'GR'
                description='Geohydrologisch model Groningen - 2008';   % Actueel (gebaseerd op REGIS II.1)
            case 'LB'
                description='Geohydrologisch model Limburg - 2008';
            case 'OV'
                description='Geohydrologisch model Overijssel - 2008';
            case 'DR'
                description='Geohydrologisch model Drenthe - 2005';     %Gedateerd (gebaseerd op REGIS II.0)
            case 'FL'
                description='Geohydrologisch model Flevoland - 2005';
            case 'FR'
                description='Geohydrologisch model Friesland - 2005';
            case 'GE'
                description='Geohydrologisch model Gelderland - 2005';
            case 'GR'
                description='Geohydrologisch model Groningen - 2005';
            case 'NB'
                description='Geohydrologisch model Noord-Brabant - 2005';
            case 'NH'
                description='Geohydrologisch model Noord-Holland - 2005';
            case 'OV'
                description='Geohydrologisch model Overijssel - 2005';
            case 'UT'
                description='Geohydrologisch model Utrecht - 2005';
            case 'ZL'
                description='Geohydrologisch model Zeeland - 2005';
            case 'ZH'
                description='Geohydrologisch model Zuid-Holland - 2005';
        end
end

%% difference between UNIX and WINDOWS
if ispc
    sx=sprintf('%.2f;',xRD); sx=sx(1:end-1);
    sy=sprintf('%.2f;',yRD); sy=sy(1:end-1);
    ids=';';
else
    sx=sprintf('%.2f\\;',xRD); sx=sx(1:end-2);
    sy=sprintf('%.2f\\;',yRD); sy=sy(1:end-2);
    ids='\;';
end

%% Construct URL
s=...
{'http://www.dinoloket.nl/dinoLks/map/SectionImage', '';
'?xCrds=',  sx;
'&yCrds=',  sy;
'&ids=',   ids;
'&dataType=', 'hge';
'&modelName=',   modelName;
'&description=', description;
'&gridResolution=' '100';
'&endDepth=', ''};

% replace spaces by '%20'
for i=1:size(s,1)
    J=findstr(' ',s{i,2});
    if ~isempty(J)
        for j=length(J):-1:1
            s{i,2}=[s{i,2}(1:J(j)-1) '%20' s{i,2}(J(j)+1:end)];
        end
    end
end

s=s';

%% Send URL
URL=[s{:}];

%web(URL,'-browser');

% I would like to retrieve the cross section direclty into a Matlab image
% however, without more details about how to request it, it seems
% difficult. Thanks to the documentation there is no such difficulty with
% respect to Google Maps figures, where the picture format is included in the request.
%  =====  [status,f]=urlwrite(URL, 'samples.html');
% In that case we can retrieve the figure
A    = imread( URL,'png');               % get it directly into matlab
info = imfinfo(URL,'png');               % get associated image info
image(xLim,yLim([2 1]),A); colormap(info.Colormap); % Plot it using the colormap from the info
% Here we can save the website. (right button, save as). We then obtain all
% files consituting the site, which includes the png files of the cross
% section and separately of its legend.
% TO 110514

%% alternative, Matlab's own browser (not so good)
%web([s{:}]);

end

function name=pclean(name)
if ischar(name)
    [P N E]=fileparts(name);
    if isempty(E),
        name=[name '.kml'];
    end
    if ~exist(name,'file'),
        name=0;
    end 
else
    name=0;
end
end


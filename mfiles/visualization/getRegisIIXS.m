function URL=getRegisIIXS(kmlpathfile,im)
%GETREGISIIXS Retrieves a cross section from the TNO dino loket website
%
% ToDo: This file has some trouble with the shell needs adaptation
%
% USAGE:
%   URL=getREGISXS(kmpathfile)
%
% REGIS is the national groundwater database of the Netherlands
%   this function gets a cross section anywhere through the Netherlands
%   directly from this database into your browser. The input file kmlfile
%   is a pathfile you get from GooogleEarth. This pathfile is scanned, the
%   coordinates are extracted and translated to the Dutch National Coordinate
%   system and inserted into the URL, after which the cross section is
%   requested from the REGIS databse. You can save the obtained page in your
%   browser yourself directly from your browser. This yields a html file that
%   will regenerate the cross section at any moment. Of course, you can also
%   make a screendumpt if desired or a selection for your report.
%
% See also: getRegis getDinoloketXSec
%
% TO 100612 110515 120409

[xRD,yRD]=kmlpath2rd(kmlpathfile);

modelName={
    'gws';
    'gws_li_geo';
    'gws_li';
    'HGM-LI-HYP'};

if nargin>1
    im=min(length(modelName),max(1,im));
else
    im=1;
end

url={
     'http://regisloket.nitg.tno.nl/dino_section/sectionServlet?';
     'gridResolution=100';
     '&xCrds=';
     '';
     '&yCrds=';
     '';
     '&modelName=';
     ''
     };

%% difference between UNIX and WINDOWS
if ispc
    sx=sprintf('%.2f;',xRD); sx=sx(1:end-1);
    sy=sprintf('%.2f;',yRD); sy=sy(1:end-1);
    sx=sx(1:end-1);
    sy=sy(1:end-1);
else
    sx=sprintf('%.2f;',xRD); sx=sx(1:end-1);
    sy=sprintf('%.2f;',yRD); sy=sy(1:end-1);
end

url{4}=sx;
url{6}=sy;
url{8}=modelName{im};

URL=[url{:}];

web(URL,'-browser');

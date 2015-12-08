function [A,xRD,yRD]=getRegisXS(varargin)
%GETREGISXS: get a cross geo(hydro)logical cross section through the Netherlands (DEPRECATED in favor of GETDINOXSEC)
%
% [A,xRD,yRD]=getRegisXS(kmlfilename)
% [A,xRD,yRD]=getRegisXS(E,N)
%
% Get REegis cross section along lines of coordinates in E,N (from path Google
% Earth either given as wgs coordinates or as kml file)
% A is the image
% xRD,yRD coordinates in Dutch National System
%
% Frans Schaars, Theo Olsthoorn 091215 dprecated as off: 1104027

if nargin==1
    if any(varargin{1}=='*')
        d=dir(varargin{1});
        for i=1:length(d)
            fprintf('Downloading Regis XS for file <<%s>>\n',d(i).name);
            [A,xRD,yRD]=getRegisXS(d(i).name);
        end
        return
    else
    [P,fname]=fileparts(varargin{1});
    [xRD,yRD]=kmlpath([fname '.kml']);
    end
elseif nargin==2
    E=varargin{1}; N=varargin{2};
    [xRD,yRD]=wgs2rd(E,N);
else
    getRegisXS
end

xCrds=sprintf('%.0f;',xRD);
yCrds=sprintf('%.0f;',yRD);


URL=['http://regisloket.nitg.tno.nl/dino_section/sectionServlet?gridResolution=100',...
    '&xCrds=', xCrds, ...  % for example 104736;111056;115594;121699;120997;129532
    '&yCrds=', yCrds, ...  % for example 490689;485233;487286;485611;494525;494147
    '&modelName=gws&gridResolution=100'];


A=imread(URL,'png');
figure; image(A);  set(gca,'visible','off');


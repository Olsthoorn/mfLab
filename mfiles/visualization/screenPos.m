function pos = screenPos(fh,fv)
%SCREENPOS gives position of figure, so that fig is centralized and its size is given fraction of screen size
%
% Example:
%    pos = screenPos(fracHor,fracVert)
%    figure('position',screenPos(0.85,0.15));
%    figure('position',screenPos(0.75));
%
% TO 130328

if nargin<2, fv=fh; end

if ~isnumeric(fh) || fh<=0.01 || fh>1 || ~isnumeric(fv) || fv<=0.01 || fv>1
    error('%s: argument(s) (screen fraction) must be between 0.01 and 1',mfilename);
end

p   = get(0,'screenSize');

pos = round([1+(1-fh)*p(3)/2, 1+(1-fv)*p(4)/2  fh*p(3) fv*p(4)]);


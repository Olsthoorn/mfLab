function h = circle(x0,y0,r,varargin)
%CIRCLE plots a circle
%
%  Example:
%     h = circle(x0,y0,r,varargin)        % draws circle
%     h = circle(x0,y0,[rx,ry],varargin); % draws ellipsoid
%
% TO 121101

n = 2*36;
x = x0+ r(  1) * cos((0:n)*2*pi/n);
y = y0+ r(end) * sin((0:n)*2*pi/n);

if nargin<3, varargin = {'visible','on'}; end

if ~strmatchi('facecolor',varargin) && ~strmatchi('edgecolor',varargin)
    h = plot(x,y,varargin{:});
else
    h = patch(x,y,varargin{:});
end


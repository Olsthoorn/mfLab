function aY = aYear(varargin)
%%AYEAR length of an avg year: =(datenum(20001,1,1)-datenum(-4000,1,1))/6000
%
% USAGE: aY = aYear(span)
%  span can be [y1 y2]
%  or numberOfYears
%
%  computes the average year length between y1 and y2
%  or between this year and +/- numberOfYears/2
%
% example:
%   aYear
%   aYear(500)
%   aYear(1780,2020)
%   aYear( [1780 2010] )
%   aYear(1780:2020)
%   aYear(2020:1780)
%   aYear(-500:1500)

% TO 131220

if nargin<1
        aY= 365.2425;
else
    if nargin==2,
        span = [varargin{1} varargin{2}];
    else
        span = varargin{1};
    end
    
    if numel(span)<2
        y = datevec(now); y=y(1);
        span = y+round([-span +span]/2);
    else    
        span = sort(round(span));
        span = span([1 end]);
    end
    aY = (datenum(span(2),1,1)-datenum(span(1),1,1))/diff(span);
end
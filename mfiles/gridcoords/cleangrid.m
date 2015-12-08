function xGr=cleangrid(xGr,dxmin)
%CLEANGRID remove columns (or rows) smaller than given value from grid
%
% Example:
%   xGr=cleangrid(xGr,dmin); %   xGr Gridline coordinates
%   yGr=cleangrid(yGr,dmin); %   yGr Gridline coordinates
%
% only used in mf_adaptPipe
%
% See also: makegrid
%
%   TO 100327

if nargin<2, error('not enough input arguments in cleangrid'); end

xGr=unique(xGr);

N=length(xGr); dx=abs(diff(xGr));

i=1;
while i<length(dx);
    if i<length(dx) && dx(i)<dxmin-eps
        dx(i+1)=dx(i)+dx(i+1);
        dx(i)=[];
    else
        i=i+1;
    end
end
xGr=xGr(1)+cumsum(dx);

fprintf('Number of  cells reduced from %d to %d by removing cells smaller than %g from vector.\n',N-1,length(xGr)-1,dxmin);


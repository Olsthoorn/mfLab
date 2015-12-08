function L = getTop(o,zoneArray,zone)
%GETTOP gets top of where zoneArray==zone
%
% Usage:
%    L = gr.getTop(zoneArray,zone)
% INPUT:
%    zoneArray of gr.size. May be a logical array
%    zone one or more zone numbmers or true (default)
% OUTPUT:
%    logial array of gr.size with one non-zero line telling which cells are
%    the top of zone.
%
% Example:
%    gr = gridObj(0:1:50,[-0.5 0.5],0:-1:-20);
%    zoneArray = gr.ZMlay <= bsxfun(@times,ones(1,1,gr.Nz), 10*sin(gr.xm/diff(gr.xm([1 end]))/5 * (2*pi)));
%    figure; spy(zoneArray)
%
%See also: bcnZone
%
% TO 120401

if ~all(o.size==size(zoneArray))
    error('%s: size of grid [%d %d %d] must match size of zoneArray [%d %d %d]',...
        o.size,size(zoneArray));
end
if nargin<3
    if  logical(zoneArray)
        zone=true;
    else
        zone=true;
        zoneArray = zoneArray~=0;
    end
end

L = ismember(zoneArray,zone);
L = diff(cat(3,false(size(L(:,:,1))),L),1,3) > 0;

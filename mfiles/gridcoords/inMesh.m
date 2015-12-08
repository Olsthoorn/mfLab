function [jc,ic,zmR,ds]=inMesh(x,y,xP,yP,zP)
%INMESH Puts polyline into a mesh, and yields mesh indics [jc=rows,ic=cols].
%
% Example:
%   [jc,ic[,zmR,ds]]=inMesh(x,y,xP[,yP[,zP]])
%
% OUTPUT:
%   [jc,ic]  the list of index pairs of cells intersected by contour xP,yP
%   zmR      z-value of contour in intersected cell (interpolated along contour)
%   ds       length of cell intersection
%   x,y      mesh grid coordinates
%   xP,yP    polyline coordinates
%   zP       is contour elevation at xP and yP
%
%   This function is suitable to generate rivers along arbitrary tracks
%   in the MODFLOW context
%
%  TODO: Upgrade this function to work with the gridObj TO130428
%
% See also:  inpoly inpolygon inpolyz hit between fallsIn above below after beyond outside
%
%   TO 070703  070706


% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin==3 && size(xP,2)>=2
    if size(xP,2)>2
        zP=xP(:,3);
    end
    yP=xP(:,2);
    xP=xP(:,1);
end
if ~exist('zP','var'),
    zP=zeros(size(xP));
elseif length(zP)~=length(xP)
    zP=ones(size(xP))*zP(1);
end

[xR,yR]=getinterface(x(:)',y(:),1,xP,yP);
ds=sqrt(diff(xR).^2+diff(yR).^2); s=cumsum(ds);

xmR=0.5*(xR(1:end-1)+xR(2:end));
ymR=0.5*(yR(1:end-1)+yR(2:end));
ic=zeros(size(xmR));  %ColNrs
jc=zeros(size(ymR));  %RowNrs
for i=1:length(xmR)
    ic(i)=find(x<xmR(i),1,'last');
    jc(i)=find(y<ymR(i),1,'last');
end

zmR=s/s(end)*(zP(end)-zP(1))+zP(1);


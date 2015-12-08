function setmulticolormap(clrstr)
% SETMULTICOLORMAP sets colormap and cmap so that multiple color schemes
%   can be used on the same axes.
%
% USAGE
%   setmulticolormap(clrstr)
%     clrstr=struct array with following fields
%         ax          % axis handel
%         range       % data value range
%         map         % the colomap
%
%   Even though each axis hat its own clim, they all share the same colormap
%   because that is a figure property. So we have to manipulate the
%   colormap to make sure each axis plots int the right part of it so that
%   individual color patterns may be used.
%
%   sx = array axes for on which colors are to be set
%   range = vectorN times 2
% TO 110423  040502

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

cL=zeros(length(clrstr),5);

for i=1:length(clrstr)
    cL(i,:)= [clrstr(i).range(1) clrstr(i).range(2) NaN NaN size(clrstr(i).map,1)];
    sz=size(clrstr(i).map);
    if sz(2)>3,
        clrstr(i).map=clrstr(i).map';
    end
end

cL(:,4)=cumsum(cL(:,end));
cL(:,3)=cumsum(cL(:,end))-cL(:,end)+1;
L=sum(cL(:,end));

cmap=NaN(L,3);
for i=1:length(clrstr)
   cmap(cL(i,3):cL(i,4),:) = clrstr(i).map;
   set(clrstr(i).ax,'clim', mf_clim(cL(i,1),cL(i,2),cL(i,3),cL(i,4),L), 'nextplot','add');
   %if i>0, set(clrstr(i).ax,'color','none'); end
   set(clrstr(i).ax,'color','none');
end

colormap(cmap);

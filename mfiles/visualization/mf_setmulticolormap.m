function mf_setmulticolormap(clrObj,varargin)
%MF_SETMULTICOLORMAP sets colormap and cmap so that multiple color schemes can be used on the same axes.
%
% USAGE
%   mf_setmulticolormap(clrObj,varargin)
%     clrObj=struct array with following fields
%         ax          % axis handel
%         range       % data value range
%         map         % the colomap
%         varargin     % addtional axis options in name-value pairs
%
%   Even though each axis hat its own clim, they all share the same colormap
%   because that is a figure property. So we have to manipulate the
%   colormap to make sure each axis plots int the right part of it so that
%   individual color patterns may be used.
%
%   sx = array axes for on which colors are to be set
%   range = vectorN times 2
%
%  Seldom used, because it complicats things. Whenever possible use only
%  one colormap per figure. Use clim to limit the map to specifi values.
%  Use true RGB color for objects and images that should not be affected by
%  any colormap.
%
% TO 110423  110502 111010 

cL=zeros(length(clrObj),5);

for i=1:length(clrObj)
    cL(i,:)= [clrObj(i).range(1) clrObj(i).range(2) NaN NaN size(clrObj(i).map,1)];
    sz=size(clrObj(i).map);
    if sz(2)>3,
        clrObj(i).map=clrObj(i).map';
    end
end

cL(:,4)=cumsum(cL(:,end));
cL(:,3)=cumsum(cL(:,end))-cL(:,end)+1;
L=sum(cL(:,end));

cmap=NaN(L,3);
for i=1:length(clrObj)
   cmap(cL(i,3):cL(i,4),:) = clrObj(i).map;
   set(clrObj(i).ax,'clim', mf_clim(cL(i,1),cL(i,2),cL(i,3),cL(i,4),L), 'nextplot','add');
   set(clrObj(i).ax,'color','none');
   if iscell(varargin)
       set(clrObj(i).ax,varargin{:});  % additional options
   end
   set(clrObj(i).ax,'xlim',get(clrObj(1).ax,'xlim'));
   set(clrObj(i).ax,'ylim',get(clrObj(1).ax,'ylim'));
end

linkaxes([clrObj.ax]);

colormap(cmap);




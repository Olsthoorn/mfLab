function cmap=mf_setColormap(ax,ranges,submaplengths,mapnames)
%cmap=mf_setColormap(ax,ranges,submaplengths,mapnames)
% set figure colormap such that each axis uses the cmap with given lengththe
% for the values in ranges
% ax is an array of axis handles
% ranges is an array with values that span the color range
% lengths if the length of each color axis
% mapnames is a set of colormap names
% There is no ouput except the colormap produced. This will however already
% be set in this function.
%
% TO 110327
%

if ~iscell(mapnames), mapnames={mapnames}; end
if ~iscell(ranges),   ranges  ={ranges};   end

N=length(ax);

cL=NaN(N,5); % data for each sub color axis

if length(ranges)~=length(ax),
    error('the number of ranges(%d) in mf_setColormap must equal the number of axes(%d)',length(ranges),length(ax));
end

for i=1:length(ranges)
    iR=min(i,length(ranges));
    iL=min(i,length(submaplengths));
    cL(iR,:)=[min(ranges{iR}),max(ranges{iR}),NaN,NaN,submaplengths(iL)];
end

% position of submap in total map
cL(:,4)=cumsum(cL(:,end));
cL(:,3)=cumsum(cL(:,end))-cL(:,end)+1;
L=sum(cL(:,end)); % total length of colormap

cmap=NaN(L,3);  % this will be the colormap to fill in next

% use mflab function mf_clim to set values of clim for each submap to make
% sure that the desired colors use the right color from the map
for i=1:size(cL,1),
   iM=min(i,length(mapnames));

   cmap(cL(i,3):cL(i,4),:) = eval([ mapnames{iM} sprintf('(%d)',cL(i,end))]);
   set(ax(i),'clim', mf_clim(cL(i,1),cL(i,2),cL(i,3),cL(i,4),L), 'nextplot','add');
   if i>0, set(ax(i),'color','none'); end
end

colormap(cmap); % set the colormap

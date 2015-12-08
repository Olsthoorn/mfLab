function plotobj(obj,faceclr,edgeclr,alpha)
%PLOTOBJ plots an object given faceclr, edgeclr and object's tranparancy
%
% Example: 
%    plotgrid(obj [,faceclr [,edgeclr, [alpha]]])
%
% Used by: 
%      mflab/examples/KwelplasjesHarmen
%      mflab/examples/mf2005/AWD_gallery_16/DrainsReportTO
%
% ToDo: replace by xsConf obj
%
% TO 100115

if nargin==1, faceclr='b'; end
if nargin<=2, edgeclr=faceclr; end
if nargin<=3, alpha=0.25; end

hold on

if nargin==1, clr='c'; end

for i=1:length(obj)
    patch(obj(i).x,obj(i).y,ones(size(obj(i).x)),'facecolor',faceclr,'edgecolor',edgeclr,'facealpha',alpha);
end


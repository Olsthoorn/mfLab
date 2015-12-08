%% Dredging of the river Vecht (Netherlands, Prov. Utrecht) 2011
% This example deals with soil mechanics, i.e. the risk of bursting of soil
% by over pressure from below (see introduction).
% The second objective is to demonstrate the use of the xsConfObj. This
% object facilitates working with cross sections by allowing to define one
% as a set of zone-layer combinations and specified properties.
% See the help of xsConfObj and the workbook of the current
% example. This workbook stores the information of the cross section in two
% worksheets, "config" and "materials". The first specifies the layout of
% the cross section, the second specifies the materials used together with their
% properties and dimensions.

%% Introduction
% To remove the pollution and the nutrients that had accumulated for over a
% century, the bottom of the river Vecht (province of Utrecht, the Netherlands)
% was dredged as from 2011. However,the water level in this river is higher
% than the land and its maintained ditch water level on either side. This is
% often the case in the Netherlands due to subsidence that occured over several
% last centuries as a consequence of drainage.
% Removal of the the hydraulic resistance of the river bottom by dredging
% causes a water-pressure rise under the adjacent land. This increased
% groundwater pressure could break the soil layers with unmanageable
% groundwater upwelling as a feared consequence.
% This being unacceptable and the standing standard safety rules not allowing
% the planned dredging because of this risk, forced to model the
% process to evaluate the actual safety and to follow the progress as
% dredging proceeds as a prerequisite for the dredging project to proceed.
%
% The modeling was done with this model in mfLab.
%
% The model is a vertical cross section perpendicular to the river.
% Different cross sections along the Vecht can readily be input using the
% convenient config and material worksheets in the workbook, while allowing
% a detailed finite differnce grid for accurate meter-scale results.
%
% The modeling is done by finite elements using fdm2.m in the mfLab
% environment. MODFLOW was not used.
%
% TO 001203 JB+TO 101123 TO 101201 101222 120804

clear variables; close all

basename='VechtDredging';

%% Request a xsConfObj
Conf = xsConfObj(basename,'Config','Materials');

%% Grid

dxMin=2;  % standard cell width

xGr = Conf.xL(1):dxMin:Conf.xR(end); % uses left and right coords of zones

gr = gridObj(xGr,[-0.5 0.5],XS(Conf.array2D('Z'     ,xGr)));

%% Get necessary model arrays using material list and config
STRTHD = Conf.array3D(gr,'head'   ,'k');  % initial heads no NaNs
HK     = Conf.array3D(gr,'kh'     ,'k');  % initial horizontal conductivity
VK     = Conf.array3D(gr,'kv'     ,'k');  % initial vertical conductivity
topHd  = Conf.array3D(gr,'tophead','k'); % initial includes NaNs

VK(isnan(VK)) = 0;
HK(isnan(HK)) = 0;

%% Find which cells are wet
WET    = min(max( (STRTHD-gr.Z(:,:,2:end))./gr.DZ,0),1); % wet is when saturation > 0

%% Boundary array, determining which cells are active and which are fixed
IBOUND = Conf.array3D(gr,'IBOUND','k');  % initial IBOUND all ones (active)
IBOUND( ~WET) = 0;                    % dry cells inactive

IBOUND( mf_find(WET,'first',1) & ~isnan(XS(ones(gr.Nz,1)*topHd)) )=  -1; %fix  top wet cells unless isnan(topHd)

%% Adapt hydraulic conductivities with degree of wetting
HK     = HK.*WET;
VK     = VK.*WET;

save underneath Conf WET
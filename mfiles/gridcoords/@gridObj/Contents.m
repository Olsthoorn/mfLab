% @GRIDOBJ
%
% Files
%   BCN2Layer           - layer = gr.BCN2Layer(BCN,Column)
%   bcnZone             - gridObj/bcnZone: [BCN,PNTSRC] = gr.bcnZone(basename,type,zoneArray,zoneVals,concArray)
%   cutout              - [VAR VARtype] = gr.cutout(VAR,VARtype,Ix,Iy,Iz)
%   envHead             - STRTHD = gr.envHead(STRTHD,STCONC,drhodc,rho0);
%   fillAquitards       - h = fillAquitards(ax,gr,VKCB,iRow,varargin)
%   fillLayers2D        - h = fillLayers2D(ax,gr,K,iRow,ILay,varargin)
%   getBCN              - [bcn,varname]= gr.getBCN(o,type,var) -- 3D array with boundary conditions
%   gridObj             - gr = gridObj(xGr,yGr[,zGr [,LAYCBD[, MINDZ[, AXIAL]]])
%   hydrostaticPntwHead - [STRTHD,p] = gr.hydrostaticPntwHead(STRTHD,STCONC,drhodc,rho0);
%   plotgrid            - gr.plotgrid(ax,clr,well,figname,figcoords)
%   plotLayers          - h = o.plotLayers(ax,layers,C,varargin) -- Plots layers of grid
%   plotLocations       - gridObj/startLoc: h = gr.plotLocations(locations,varargin)
%   plotMesh            - h = o.plotMesh(ax,[color],varargin) -- Plots wireframe of grid
%   plotMeshUpdate      - h = gr.plotMeshUpdate(h,Var,varargin)
%   plotWells           - hw = gr.plotWells(ax,wells)  -- plots   wellscreens in 3D
%   plotXS              - gr.plotXS(ax,iy,well,varargin)  plot wells on X-section.
%   plotXSec            - gridObj/plotXSec: h=gr.plotXSec(jRow [,layers [,varargin]])
%   plotYS              - gr.plotYS(ax,ix,well,varargin)  plot wells on Y-section.
%   plotYSec            - plot cross section along y-axis
%   relloc2model        - modelLoc = relloc2model(o,locations)
%   setMinDZ            - gr.setMinDZ(MINDZ,layers)
%   startLoc            - gridObj/startLoc: startloc = gr.startLoc(basename,zoneArray,zoneVals)
%   well                - [well,WEL,PNTSRC]=gr.well(basename,HK,well_or_wellSheetName)
%   XSlayout            - MF_CONF: Gets grid values of the parameter parnam from specified zones

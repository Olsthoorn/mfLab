clear all;close all;clc

load basename;                % needed for B,H files and for nper.
load underneath;              % containg tje gridObj

res=resObj(basename,gr);      % setup the results object

res.cRange=[0 1];             % concentration/density-multiplier
res.lim;                      % sets x,y,z lim from gridobj;loads c-lim from cRange

res.strmRange=-32:8:32;       % these are the values for each flowline
res.strmPlotarg={};           % containing the plot vargins, all patchproperties may be used

res.equiConcRange=[0 .1 .2 .5 1];% these are the values for each flowline
res.equiConcPlotarg={};       % containing the plot vargins, all patchproperties may be used

res.concSliceX=[0];           % a concentration slice/plane on x-ax
res.concSliceY=[1 2];         % a concentration slice/plane on y-ax
res.concSliceZ=[];            % a concentration slice/plane on z-ax
res.concSlicePlotarg={'FaceAlpha',0.6,'EdgeAlpha',0};% containing the plot vargins, all patchproperties may be used

res.concfX=[];                % a concentration slice/plane on x-ax, using contourf
res.concfY=[];                % a concentration slice/plane on y-ax, using contourf
res.concfZ=[-50];             % a concentration slice/plane on z-ax, using contourf
res.concfArg=[8];             % set contourf option; n or v in the contourf helpfile, default:8
res.concfPlotarg={};          % containing the plot argins, all patchproperties may be used

for i=1:res.nper
    res.update(i);
    % you can use your getframe here. todo: a makemovie function within
    % resObj
end
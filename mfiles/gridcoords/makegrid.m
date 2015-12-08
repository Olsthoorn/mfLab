function [XGR,YGR]=makegrid(well,aroundAll,~,dmax,dmin,nstep)
%MAKEGRID geneate a grid that is refined around given wells
%
% Example:
%   [XGR,YGR,well]=makegrid(well,aroundAll,aroundEach,dmax,dmin,nstep)
%
%    First version by Ruben Calje, then simplified by TO to more radical grid generation.
%    Finally improved by Mark van der Valk for his MSc.
%
%    The basic grid has uniform cell size dmax and covers all wells plus
%    a rectangular space around the well set of size "aroundAll". 
%    Each well is surrounded by a fine  grid starting with cell size dmin
%    increasing to the cell size of the overall grid, dmax in nsteps steps.
%    These meshes are then merged; columns and rows smaller than dmin are
%    set to the average of the two closeby grid lines, eliminating one.
%    "Ugly" leftover gridlines are finally removed.
%    No coordinates are needed to generate this mesh. The variable dmin, dmax,
%    aroundAll, aroundEach and nsteps are sufficient, together with the
%    coordinates of the wells in the well array "well".
%
% RCalj?  2009
% MvdValk 2010
% TO 091117 091124 120824 150119


if nargin<2, [XGR,YGR] = selftest(); return; end

fprintf('Generating computation FDM grid for this model ....\n');

%% aroundAll  is the space of model around the bounding box of all wells
%  aroundEach is the space of refind grid around each well  
if dmax<=dmin,               error('%d: dmax<=dmin',mfilename); end
if nstep<1,                  error('%d: factor<1',  mfilename); end

%% Round off well coordinates to center of cell with size dmin
xw = round([well.x]/dmin)*dmin;
yw = round([well.y]/dmin)*dmin;

% A general coarse network
XGR = min(xw-aroundAll):dmax:max(xw+aroundAll);
YGR = min(yw-aroundAll):dmax:max(yw+aroundAll);

%% subgrid to be placed around each well

dx = logspace(log10(dmin),log10(dmax),nstep); % increasing widths
d  = cumsum([dx(1)/2 dx(2:end)]); % subgrid positive coordinates
L  = d(end); % extent of half grid

subgrid=[-d(end:-1:1)  d]; % subgrid span

%% Delete coarse grid within subgrid of the wells because that of the wells is finer
for iW=1:length(well)
    XGR=XGR(XGR<well(iW).x-L | XGR>well(iW).x+L);
    YGR=YGR(YGR<well(iW).y-L | YGR>well(iW).y+L);
end

%% Add subgrid lines to coarse grid lines
Nx = length(XGR);
Ny = length(YGR);
Ns = length(subgrid);
Nw = numel(well);

XGR = [XGR NaN(1,Nw*Ns)];
YGR = [YGR NaN(1,Nw*Ns)];
    
for iw=1:length(well)
    XGR(Nx+ (iw-1)*Ns+(1:Ns)) = well(iw).x+subgrid;
    YGR(Ny+ (iw-1)*Ns+(1:Ns)) = well(iw).y+subgrid;
end

%% Remove obsolete grid lines
XGR=cleangrid(unique(round(XGR*100)/100),dmin);
YGR=cleangrid(unique(round(YGR*100)/100),dmin);

fprintf('.... Grid construction finished.\n');

%% Clean grid (removing cells smaller than dmin)
function XGR=cleangrid(XGR,dmin)

k=0;
while 1
    Dx=diff(XGR);
    
    % which inner grid point has the smallest cell either left or right of it?
    minDx = min([Dx(1:end-1); ...  % min dx of cells left to inner grid lines
                 Dx(2:end)]);      % min dx of cells right to inner grid lines
    
    minminDx = min(min(Dx));       % smallest overall cell width
    
    % now get grid line to be removed
    if rem(k,2)==0  % not always choose same side
        imin= find(minDx == minminDx,1,'first') + 1;  % grid point
    else
        imin= find(minDx == minminDx,1,'last' ) + 1;  % grid point
    end
    
    % remove it or ready?
    if minminDx<dmin
        XGR(imin)=[];
        k=k+1;
    else
        return;
    end
end

function [XGR,YGR] = selftest()
    %% selftest -- self test for makegrid
    % generates some wells from workfile 'testMakegrid.xls', sheet 'wells'
    % which have random coordinates and then generates a grid and plots it.
    % Each run yields different locations due to rand function in excel
    % when computing coordinates of wells.
    % TO 150114

    basename   = 'BTOww';
    sheetNm    = 'Wells';
    well       = wellObj(basename,sheetNm);
    dmin       =   2;
    dmax       = 100;
    aroundAll  =  500;
    aroundEach = dmax;
    nstep      = 11;
    [XGR,YGR] = makegrid(well,aroundAll,aroundEach,dmax,dmin,nstep);
    
    hold on; xlabel('x [m]');  ylabel('y [m]'); title('test makegrid')
    plotgrid(XGR,YGR,'c',1,well)
    plot([well.x],[well.y],'ro');
 

function [Values,IDX] = interpLaplacian(o,ixy,vertexData,onIdx)
% GRIDOBJ/INTERPLAPLACIAN interpolate contour values to internal area by
% solving laplacian
%
% USAGE: [Values,inIdx] = gridObj.interpPlaplacian(xy,vertexData[,onIdx])
%
% Inputs:
%    gridObj must be an instance of a grid object.
%    xy are the coordinates of a polygon, being the circumference of an area
%    into which the contour values are to be interpolated.
%
%    if onIdx is unspecified,
%      then vertexData is an array of which each column is a data type to be
%      interpolated and each row coincides with a vertex of the area
%      boundary and the number of rows of vertexData must equal the number of
%      rows of the boundary vertices.
%    if onIdx is specified,
%      then vertexData must correspond to the onIdx cell indices, and,
%      therefore, size(vertexData,2) should equal numel(onIdx).
%    if onIdx is unspecified, onIdx and the vertexDat on cells onIdx will
%      be computed internally.
%
%    onIdx are the cell indices intersected by the polygon.
%
% Algorithm
%    The Laplacian interpolation is carried out using a 3D steady-state
%    finite difference groundwater model by setting Kx=Ky=1 and Kz=0, fixed
%    head at the bounaries equal to the vertex data columns and zero external
%    flows. Each layer in the model represents one vertex column, i.e. a
%    different data type.
%
% TO 130820

ixy = ixy(:)';
if numel(ixy)~=2
    error('size iXY must be 2 (indices of x and y columns)');
end

if nargin<2    
    error('3 inputs required: index of [x y] vertexDat and onIdx');
end
   
xy = vertexData(:,ixy);

ilay = floor(median(onIdx(:)-1)./o.Nxy)+1;  % choose layer of majority of lineObj cells

if nargin<3
    % generate a line through the model with zRel=0 (second argument)
    P = o.lineObjects(xy,0);
    onIdx = [P.Idx];

    % make sure that onIdx refers to layer 1
    onIdx = onIdx - o.Nxy * floor((onIdx-1)/o.Nxy);

    %% Interpolate all values to their grid locations

    % cumulative line, xy, length at vertices
    sL = [0; cumsum(sqrt(sum(diff(xy,1,1).^2,2)))];

    % cumulative length along the intersected cells
    % to the mid of the intersection with the cell
    sCell = [0 cumsum([P_.L])];
    sCell = 0.5*(sCell(1:end-1)+sCell(2:end))';
    % elevation is known to P don't need to compute it

    % Interpolate vertices avlues to intersected cell centers,
    % keep worksheet order
    vertexData = interp1(sL,vertexData,sCell);
 
else

    if ~isvector(onIdx), error('%s: onIdx must be a numerical vector of indices',mfilename); end

    if numel(onIdx) ~= size(vertexData,1)
        error('%s: numel(onIdx) must equal size(vertexData,1) of onIdx is specified',mfilename);
    end
    
    % make sure that onIdx refers to layer 1
    onIdx = onIdx - o.Nxy * floor((onIdx-1)/o.Nxy);
end

[In,On] = inpolygon( o.Xm,o.Ym, xy(:,1),xy(:,2));

inIdx = find(In| On);

Ncol = size(vertexData,2);

IB     = zeros(o.Ny,o.Nx); IB(inIdx)=1; IB(onIdx)=-1;                
IBOUND = bsxfun(@times,IB,ones(1,1,Ncol));

Phi    = zeros(o.Ny,o.Nx,Ncol);
FQ     = zeros(o.Ny,o.Nx,Ncol);
IDX    = NaN(numel(inIdx),Ncol);

for iCol=Ncol:-1:1
    offset = (iCol-1)*o.Nxy;
    Phi(onIdx+offset) = vertexData(:,iCol);
    IDX(:,iCol) = inIdx + offset;
end

if Ncol==1
    Phi = fdm2(o.xGr,o.yGr,    1,1,  IBOUND,Phi,FQ);
else
    zGr = 0:-1:-Ncol;
    Phi = fdm3(gridObj(o.xGr,o.yGr,zGr),1,1,0,IBOUND,Phi,FQ);   % arguments:  gr,kx,ky,kz,IBOUND,IH,FQ
end

Values = Phi(IDX);

IDX = IDX(:,1); % only the first column are the actual indices

% replac x and y with the values of the cells
Values(:,ixy) = [o.XM(IDX) o.YM(IDX)];

IDX = IDX + (ilay-1)*o.Nxy;

%% It is possible to compute this for all vertex data in one run
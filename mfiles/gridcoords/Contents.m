% GRIDCOORDS
%
% Files
%   above              - gets indices of zGr or zm before point a, where zGr is ascending
%   after              - Indices of cells after a where xGr or xm is assumed to ascend
%   bcnZone            - yields list data for stresses in MODFLOW and MT3DMS/SEAWAT
%   before             - gets indices of xGr or xm before point a, where xGr is ascending
%   below              - gets indices of zGr or zm before point a, where zGr/zm are ascending
%   between            - gets indices indices of xGr or xm between find(xGr>min(a,b) & xGr<max(a,b)).
%   beyond             - get indices of cells xGr or xm beyond a where xGr}xm is assumed ascending
%   cellIndex          - get glocal index of 2D or 3D array using individual coordinate indices
%   cellIndices        - get individual axes indices given global one, chosen order
%   cleangrid          - remove columns (or rows) smaller than given value from grid
%   convertLAYCBD      - converts layers from LAYCBDold to LAYCBDnew (with selftest)
%   cutoutArray        - cutout 3D subarray using indices Ix,Iy and definition of confing beds of input and output
%   cutoutBCN          - cuts out boundary WEL, GHB,CHD,DRN,RIV, ... to match coordinate indices Ix and Iy
%   cutoutZTA          - cuts out a piece of ZTA accoring to input rows and columns
%   deg2utm            - convert Lat Lon to x,y, utmzone
%   drainObj           - DRAIN class def for hydObj
%   factorspace        - get a power series from 0 to L with factor fac and n points
%   fallsIn            - gets indices indices of xGr or xm between find(xGr>min(a,b) & xGr<max(a,b)).
%   geo2grid           - maps geology given by zGeo and kGeo to grid given by zGr
%   getDinoXSec        - retrieve geo(hydro)logcal cross section from www.dinoloket.nl
%   getDinoXSecA       - retrieve geo(hydro)logcal cross section from www.dinoloket.nl
%   getInitialSalinity - generate sigmoid vertical salinity distribution
%   greatCircle        - gets distance given Lat Long for two points
%   gridsTransfer      - % convets valuesFr corresponding to source grid (xGrFr) to valuesTo accoring to xGrTo
%   gridsurf           - plots a surface using grid coordinates for x and y and center of cell values for ZM
%   gw2rd              - convert local Amsterdam Water Supply Dune Area coords to Dutch national coordinates.
%   hit                - get index of cells between a and b or if one input argument the cell in which a resides
%   IdxMatlab2Modflow  - converts global Matlab index to global MODFLOW index
%   inMesh             - Puts polyline into a mesh, and yields mesh indics [jc=rows,ic=cols].
%   inpoly             - test if poitns are in polygon
%   inpolyz            - returns logical arrray for zx plane telling which cells are in vertical polygon
%   isAquifer          - logical vector telling which model layer is an aqufier and which is not
%   itype              - for SSM package of MT3DMS and SEAWAT
%   JoinBCN            - joins a MODFLOW or MT3DMS/SEAWAT boundary condition list L of form [iPer iLay iRow iCol rest]
%   JoinLayers         - joins layer array OldLayer according to JoinArray
%   kmlpath            - gives wgs-coodinates of GoogleEarth path embedded in kml file
%   kmlpath2rd         - give xRD and yRD coordinates of GE path in kml file
%   kmlplacemarks      - gives wgs-coodinates of all placemarks in GoogleEarth path embedded in kml file
%   laycbd_check       - checks consistency if LAYCBD with number of layers in MODFLOW model
%   layer2aquif        - converts layer number (LAY+CBD) to aquifer number (LAY) using full LAYCBD
%   lineGrid           - get info on all line pieces of polyline pline intersecting a 1, 2 or 3D mesh
%   linegridObj        - aribrary line through a finite difference grid or mesh
%   makegrid           - geneate a grid that is refined around given wells
%   mergeModel         - convert Model into output Model with newCoords along dim
%   mf_cleandem        - removes leading and trailing rows and columns with NaNs from dem
%   mf_dem2grid        - make a dem using grid coordinates
%   mf_demCoarse       - generates a coarser dem (digital elevation model)
%   mf_find            - finds start or end of a zone in a multidimensional array along dimension dim
%   mf_getdemfromtiff  - read a DEM from a TIFF file
%   mf_kmlpath         - get wgs coords of GoogleEarth path embedded in kml file
%   mf_observe         - gets data fro given obvervation points
%   mf_setHFB          - Sets Horizontal Flow Barrier
%   mf_setmnwells      - puts wells in the model grid.
%   mf_setwells        - puts wells in the grid when they are specifie in the sheet with 
%   mf_Z_extend        - extends x and y of Z to nodes instead of cell centers
%   mf_zone            - gets grid values from conf and mat worksheets
%   movie3D            - generates a 3D movie of the ATES simulation
%   movie_xsec         - Makes movie of vertical cross section of model using output H
%   movie_xsec_test    - generate movie for vertical cross section of model using output H
%   oneZeta            - extracts zeta planes (fresh-salt interfaces) from struct ZETA (SWI, salt water intrusion package)
%   outside            - Get indices if xm outside a given range [a b];
%   plotConf           - plots configuration specified in sheets Config and Materials
%   plotgrid           - plots the grid lines in color clr given the coordinates xGr yGr
%   plotobj            - plots an object given faceclr, edgeclr and object's tranparancy
%   point2line         - puts point xp yp on the line given by end points of X(1 end) Y(1 end)
%   rd2gw              - converts Duth national coords to local coords of Amsterdam Water Supply Dune Area
%   rd2wgs             - converts Dutch national rd-coords (x,y) to GoogleEarth wgs coords (lon, lat) + selftest
%   readIDF            - reads IMOD's IDF file (Peter Vermeulen, Deltares)
%   RefineBCN          - refines stresses lists WEL, DRN, RIV, CHD according to SplitArray
%   RefineGrid         - refines rows or cols of a FDM model array according to splitArray
%   removeCBD          - turn all confining beds of model into model layers
%   rotate             - rotates coordinates over alfa (counter clockwise) degrees around x0,y0
%   sinespace          - generates a nice spacing based on end points
%   subArea            - cuts subarea from area given LONLIM LATLIM of area and lonlim latlim of subArea
%   utm2wgs            - converts the vectors of UTM coordinates into Lat/Lon vectors.
%   wgs2rd             - computes Dutch national coords (xRD,yRD) from GoogleEarth WGS coords (E,N)
%   wgs2xy             - adds x,y coordinates to obs struct with fields E,N usign e0,n0 as center
%   xy2wgs             - adds x,y cordinates to obs struct
%   xyzindex           - computes the cell indices of points in 1D, 2D or 3D grid
%   xyzindexTest       - tests function xyzindex

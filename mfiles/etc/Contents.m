% ETC
%
% Files
%   adaptSetPoints   - adapt set points for ATES systems
%   addcopyright     - adds or updates mfLab copyright notice
%   ADVmethod        - class definition for advection methods
%   arrayCheck       - verfies that there are no NaN's or Infs in the array.
%   avi_compress     - compresses avifile on Mac systems
%   binRange         - gets a range min(x):dx:max(x) for the data x
%   cell2list        - make list from a cell array of lists
%   cleanBinary      - removes record length info from binary files
%   colorCell        - fills cells with color
%   cols             - returns number of columns of array
%   ContourRange     - computes a suitable range of values for contouring.
%   copyright        - Copyright 2009 Theo Olsthoorn
%   datenum2excel    - converts matlab datenum to excel datenum
%   dbfread          - reads data from dbf file
%   deltaValues      - substract starting values from data in struct which must have field 'values'
%   digitiz          - digitize on screen axis using mouse and show what has been selected
%   done             - write done with of without the toc value
%   dtopBalance      - Processing water balance for spreadsheet
%   environmentalHd  - computes environmental head (Seawat boundary conditions)
%   excel2datenum    - converts excel datenum to matlab datenum
%   fileGrep         - find expr in file using regexp(fileContents,expr)
%   fprintfs         - fprintf but for multiple strings given in cell array
%   getDem           - gets dem from USGS tiffile specific for MS reservoir
%   getGXG           - computes average lowest, spring and highest groundwater level
%   getinterface     - Yields coordinates of an interface
%   getNext          - get prop value from varargin, if not present use default
%   getProp          - get property value from varargin, if not present use default
%   getTable         - reads a table into a struct with given regCols and the rest in field UserData
%   getWord          - looks of word is one of the values in varargin
%   GRFgenerator     - Gaussian Random field generator according to Hoo (2004)
%   ierfc            - repeated integral of complementary error function
%   ImageModel       - Access to properties of an image relevant to its display.
%   isaxis           - verify that var is an axis
%   ismynan          - returns true of any X==NaN
%   isotherm         - determines sorption isotherm (class def)
%   jet2             - generate useful colormap similar to jet
%   maskBud          - sets B(i).term{j}(:) to 0 for values>MASK or values<MASK for values==0
%   maskH            - masks heads H(i).values using HINACT.mat saved by mf_setup
%   maskHC           - masks H(i).(fieldName) if  isstruct(H) or H(i).(fieldName({j} if isstruct(H(i).fieldName)
%   mf_checkdir      - verifies that output has been generated, i.e. the model run terminated normally.
%   mf_cleandir      - safely removes the files named in the nam files in this directory and some more
%   mf_clim          - yields the extremes of the colormap values to be used as clim([c1 c4])
%   mf_expand        - expands matrix to full matrix
%   mf_gcounter      - puts a window with text and values on figure
%   mf_logo          - plots string www.google.com/p/mflab on lower left of figure
%   mf_Psi           - computes stream function based on FLOWRIGHTFACE or FLOWFRONTFACE
%   mf_setTime       - adds simulation time to budget struct
%   mf_stat          - statistics of input X
%   mflab            - script to run mflab in any directory
%   nl               - print a new line
%   nums2str         - converts numeric values to strings using given formats
%   pick             - picks next value from list and start with 1  after reaching end
%   piezometerObj    - class definition of piezometer Objects
%   plotshp          - plots the shape
%   readADF          - potentially reads all grid datasets of the directory one by one.
%   rm               - removes files from directory (warns to use 'delete' instead of 'rm')
%   roundn           - rounds x to n didgets.
%   setHFB           - geneate the array required for HFPB processes.
%   setmulticolormap - sets colormap and cmap so that multiple color schemes
%   sorbedMass       - computes total sorbed mass in cells IDX of grid
%   sprintfs         - sprintf for strings in cell array strs
%   strmatchi        - finds matches of (cell array of) string(s) in cell of array of strings
%   SWI_getIface     - extract fresh-salt interface from a 3D MocDense input file
%   swstate          - State equation for seawater
%   variogram        - compute a variogram from the data
%   wait4model       - wait until model that runs in background has finished
%   XS               - convert 3D array to cross section along X-axis?
%   YS               - convert 3D array to cross section along Y-axis?
%   zonebudget       - compute water budget for set of zones

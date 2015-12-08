% NHI
%
% Files
%   analyze           - analyzes layers of NHI
%   analyzeREGISgrids - analyzes regis GRIDS
%   clayers           - contours all layers where VAR is the 3D-value matrix
%   getNHIASC         - reads ASCII (ESRI) datafile, select between given coordinates
%   getNHIBCN         - gets boundary conditions/stresses from files RIV GHB DRN
%   getNHIdata        - script to read the input files of the Netherlands Hydrologic Instrument (NHI)
%   getNHIfileNm      - gets NHI file name for workbook saved in mfLab
%   getNHImeta        - gets the 6 lines of meta data from any of the NHI .ASC files.
%   getNHISCD         - retrieves the WEL stresses in NHI model on a per cel bases
%   getNHIzipfileList - gets list of zip files contained in NHI on website (www.NHI.nu)
%   locate_NHI        - starts with directory with NHIzipfiles and generate a directory and ...
%   makeArray         - generates a 3D array of the size of the NHI model or part of it given xLim,yLim
%   mkNHIbib          - starts with directory with NHIzipfiles and generates a directory with unzipped files
%   NHI               - extrac NHI data from the ASCII files in mfLab/examples/NHI/NHIascii
%   NHI2AGV           - extracts AGV model from NHI datafiles
%   NHI_1             - script to extract the data from NHI ASCII files in zipfiles at www.NHI.nu
%   NHI_getziplist    - gets list of zip files contained in NHI on website
%   NHI_read          - reads NHI datafiles and save in mat file to minimize space
%   NHI_readASC       - reads NHI ASCII datafile (ESRI), select between given coordinates
%   NHI_readmeta      - reads meta data belonging to NHI ascii data files (ESRI)
%   NHI_savezips      - saves the zip files of the NHI site
%   NHI_unzips        - downloads, unzips and saves the NHI zip files URL
%   readASC           - read NHI ascii files (ESRI ascii georefferenced format)
%   readNHImeta       - reads metaDat from NHI file
%   readSCD           - reads SCD (stress) files van NHI model (drains, rivers ghb and wells)
%   setRIVdepth       - sets layer number equal to the layer in which the RIV bottom resides.
%   setshapedir       - renames files in shapefile dir so that they all have the same basenae and
%   showshp           - plots shapfile in black by reading appropriate shapefiles
%   testshp           - tests shape file
%   unpack            - script to to get variables from ModelObj array into the workspace (through eval)

% READ
%
% Files
%   mmf_readmdl         - reads a MODFLOW model
%   namSCEN             - looks to see if PCKG is in PCKGS. If not, returns 0
%   rarray              - general MODFLOW/MT3DMS array reader
%   read_mmf            - reads a MODFLOW model into Matlab (obsolete)
%   readADV             - reads basic advection transport package file
%   readBA6             - reads MODFLOW's bas6 package file
%   readBAS             - read MODFLOW's basic bas6 package file
%   readBC6             - reads MODLFOW's basic flow package file
%   readBCF             - reads MODFLOW's basic flow package file, old version
%   readBCN             - reads MODFLOW's stress files (Boundary Condition Files)
%   readBTN             - reads MT3DMS/SEAWAT basic transport package file
%   readBud             - reads MODFLOW's budget file into a Matlab structure array.
%   readBud2            - reads MODFLOW budget file into a Matlab structure array
%   readBud6            - reads compact budet file needed by MODPATH6
%   readCFP             - reads MDOFLOW-CFP input file (Conduit flow package)
%   readCHD             - reads MODFLOW's CHD boundary package file
%   readCOC             - reads output control for CFP package (mf05)
%   readCRCH            - reads conduit flow recharge package input file
%   readDat             - reads binary MODFLOW heads or drawdown output into a struct, H
%   readDat2            - reads binary MODFLOW heads or drawdown output into a struct, H
%   readDIS             - reads MODFLOW's discretization package input file
%   readDSP             - reads MT3DMS's dispersion diffusion package input file
%   readENDP            - reads end points file produced by MODPATH6
%   readEndPoints       - read simulation end points produced by MODPATH6
%   readEVT             - reads modflow's EVT package input file
%   readHFB             - reads MODFLOW's horizontal flow barrier package
%   readHUF             - reads Horizontal Unite Flow  package file
%   readLPF             - reads layer property flow package file
%   readLVDA            - reads Layer Variabl Direction Anisotropy package file
%   readMDL             - same as readMFlab
%   readMFLAB           - reads a MODFLOW model into mfLab
%   readMFLOW           - reads an entire model into Matlab
%   readMOCM            - function mocm = readMOCM(fname,pth)
%   readMOCMAIN         - reads mocDense main file
%   readMT3D            - reads unformatted MT3DMS concentration file into struct 
%   readMT3D2           - reads unformatted MT3DMS concentration file (MT3D00n.UCN) into struct 
%   readMULT            - reads MODFLOW's multiplier file package
%   readNAM             - reads the name file
%   readOBS             - reads formatted conc obsrvation point outut 9(MT3D, Seawat)
%   readPath            - reads path line file produced by MODPATH6
%   readRCH             - reads MODFLOW's recharge package input file (RCH)
%   readshp             - reads an ERSI shape file
%   readSSM             - reads MT3DMS's source sink mixing package input file (SSM)
%   readSWI             - reads Salt Water Intrusion input file
%   readTsrPoints       - reads time series file produced by MODFPATH6
%   readVDF             - reads SEAWAT's VDF package
%   readzeta            - reads interfaces for Salt Water Intrusion package (SWI)
%   skipmodflowcomments - skips matlab comment lines, i.e. lines that start with #
%   sp_timeObj          - class def of timeObj to get time for every stress period and time step

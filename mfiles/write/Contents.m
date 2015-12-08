% WRITE
%
% Files
%   BCN_Curb       - cuts out and renumber BCN when cutting out peice of larger model
%   BCN_Struct     - makes a list of form [iPer iLay iRow iCol rest]
%   cprintf        - displays styled formatted text in the Command Window
%   fmtstr         - writes ntimes a format (i.e. %2d or so) to a string for use in sprintf
%   ftfmt2ml       - converts fortrans format to matlab format for printing
%   getExcelData   - reads info from excel worksheet
%   getLayers      - reads layer information from LAY worksheet, expands where necessary
%   getPeriods     - reads stress period information from workhsheet and expands when necessary
%   mf_run         - sctipt that runs mf_setup, mf_run is considered a more logical name than mf_setup
%   mf_setup       - generates all input files for models (Backbone of mfLab)
%   mlfmt2ft       - converts a matlab number format to a fortran equivalent.
%   setExecutables - script that sets mfLab's executables in mflab/bin
%   warray         - writes MODFLOW and MT3DMS arrays to modflow input files (general array writer)
%   writeADV       - writes basic advection transport package file
%   writeASP       - writes input file for Doherty's MODFLOW-ASP program
%   writeBAS6      - writes input file for MODFLOW's basic (BAS6) package
%   writeBCF       - writes input fiel for MODFLOW's basic flow package file (is BCF6)
%   writeBCN       - writes input file for MODFLOW's stress packages (WEL, DRN, RIV, GHB, CHD ...)
%   writeBTN       - writes input file for MT3MDS's basic transport package.
%   writeCFP       - writes the Conduit Flow Package input file
%   writeCOC       - writes intput for output control pacakge of CFP (Conduit FLow Package) version of MODFLOW
%   writeCRCH      - writes input file for conduit flow recharge package
%   writeDE4       - writes input file for MODFLOW's direct sovler package DE4
%   writeDIS       - writes input file for MODFLOW's dicretization package (DIS)
%   writeDSP       - writes input file for MODFLOW's dispersion diffusion package
%   writeETS       - writes input file for MODFLOW's ETS package
%   writeEVT       - writes input file for MODFLOW's EVT package
%   writeGCG       - writes inptu for MT3DMS's GCG generalized conjugate gradient sovler package.
%   writeHFB       - writes the input file for the HFB6 package (mf2k)
%   writeHUF       - writes input file for MODFLOW's HUF package
%   writeLAK       - writes input file for MODFLOW's lake package
%   writeLMT       - writes input file for MT3MDS LMT package
%   writeLPF       - writes input file for MODFLOW's LPF package (Layer Property Flow)
%   writeMNW       - writes input file for MODFLOW's NMW2 package
%   writeMNW1      - writs input file for multinode well (version 1) package
%   writeMPBAS     - writes input file for MODPATH6's basic package
%   writeMPTH6     - writes MOPATH6's simulation input file
%   writeMUL       - writes input file for MODFLOW's multiplier array package
%   writeNAM       - writes name files for MODFLOW, MT3DMS and SEAWAT
%   writeNWT       - writes input for MODFLOW-NWT input package
%   writeOC        - writes input file for MODFLOW's output control file (OC)
%   writeOCwords   - writes output control input file using WORDS mode
%   writePCG       - writes input for MODLFOW's PCG solver package
%   writePES       - writes iput file for PEST program (calibration for any model)
%   writeRCH       - writes input file for MODFLOW's recharge (RCH) package
%   writeRCT       - writes input file for MT3DMS's reaction package (RCT)
%   writeSEN       - writes sensitivity file MODFLOW 2000
%   writeSIP       - writes intput for MODFLOW's SIP solver package
%   writeSOR       - writes input file for MODFLOW's SOR solver package.
%   writeSSM       - writes input file for MT3DMS's sink-source mixing package SSM package
%   writeSWI       - writes input for SWI package (seawater intrusion package)
%   writeUPW       - writes input file for MODFLOW-NWT UPW pacakge (upwind flow)
%   writeVDF       - writes input file for the variable density flow package file (SEAWAT)
%   writeZON       - writes input for zone file

% Example see USGS Modflow 2000 manual, Open-File Report 00-92
%
% This is the same examples as in directory ../ex1. So look for details
% and comments regarding the basic input and methods used in this
% mf_adapt.m in the file mf_adapt.m in that directory. In this directory
% wel weil foces on calibration using PEST. Therefore, we will now only
% comment the lines specific to calibration.
%
% We will use the LPF package this time as is done in the example in the
% MODFLOW 2000 manual on page 102ff. LPF requires HK instead of TRAN and
% VKCB instead of VCONT. VK is set to a large value to minimize its
% influence on the vertical resistance between layers beyong that of VKCB.
% Instead of using the mf2k parameter functionality, we just use the
% parameters directly in mf_adapt to generate the model input files. Within
% the matlab environment, such facility is not necessary, unless one wants
% to use the built-in calibration fucntionality of mf2k (which is no longer
% available in mf2005, by the way.
%
% For instructions on how to generate the files for PEST (template files,
% instruction files and the PEST command file, see comments in
% genpestfiles.m
%
% Here we comment only how to generate input files using parameters.
%
% To use past in the mfLab context, we just need an external file with the
% parameter values. This file has the initial ones when we start and is
% replaced b PEST with the current parameter values. The file is just a
% list of the parameter values.
% 
% In mf_adapt.m we read in thes values, attribute them to parameters and
% use those parameters to generate the model arrays and boundary
% conditions. This is straightforward in matlab as everything in mfLab is
% parameterized anyway.
%
% We also specify the observation points and compute their locatins in the
% model grid, so we can readily extract the computed observations at the
% end of the run. This is done in mf_analyze. mf_analyze also writes these
% computed observations to the output file to be used by PEST.
% PEST generates a new parametr value file, invokes run.bat and so on. 
%
% To easily generate the template, instruction and control file for PEST
% see instructions in mfile genpestfiles.m.
%
% TO 100728

clear valriables; close all
basename='ex1pest';

% The parameter AFTERMFSETUP is useful for various purposes in the context
% of running PEST. It's value is the mfile to be run after mf_setup is
% ready and therefore, the models have run and produced their output.
% Hence, if this parameter is set (exists), then that mfile is
% automatically run, if not, matlab stops and waits for new commands.
% Note that mf_setup.m knows about this parameter, so the exact spelling is
% essential. If not used, just comment it out.
AFTERMFSETUP='mf_analyze';
% Te possibility of testing in mf_analyze.m of this parameter exists allows
% the logic to be different in the case we run pest of when we just run the
% model without any calibration. When we run pest we just care for the
% computed observations, and when not, we like to have contours, graphs and
% other outcomes.

%% Reading in the parameter values
try
    load(par);
catch ME
    fprintf('%s\n',ME.message);
    par.phk1  =1;
    par.phk2  =1;
    par.phk3  =1;
    par.pvkcb1=1;
    par.pvkcb2=1;
    par.pdrn  =1;
    par.prch1 =1;
    par.prch2 =1;
    save par
end

% Observation points, given in x,y coordinates and layer number
xObs = 1.0e+04 * [1.5806 2.0645 4.2097 6.4355 4.2742 2.4516 4.8226 6.2258];
yObs = 1.0e+04 * [2.1637 2.9620 2.4298 6.1345 4.6813 4.8860 3.1462 3.3509];
izObs=           [ 3     2      2      1      1      2      1      3     ];

%% Grid
NLAY=3;

xGr = 0:5000:75000;  % feet
yGr = 0:5000:75000;  % feet
zGr = [200 -150 -200 -250 -350 -450];  % Plane elevation vector

gr = gridObj(xGr,yGr,zGr);

%% Get the cells of the observation points
[ixObs,iyObs]=xyzindex(xObs(:),yObs(:),xGr(:),yGr(:));

% and use global cell index to readily colate them inthe 3D arrays upon output
IObs=cellIndex(ixObs(:),iyObs(:),izObs(:),gr);

%% Then continue generating the model arrays as usual, but now using the
% parameters 

%% IBOUND array
IBOUND = gr.const(1);
IBOUND(:,1,1:2)=-1;             % then change left column layers 1 & 2 into fixed heads

STRTHD = gr.const(0);

HK     = gr.const([par.phk1*1e-3; par.phk2*1e-4; par.phk3*2e-4]); % ft^2/c
VK     = gr.const(1.0);  % [ft/s] a high value for all vertical aquifer conductivities
VKCB   = gr.const( [par.pvkcb1*2e-6; par.pvkcb2*5e-7] );             % [ft/s] vertical K of confining unit below layer 2

%% RCH (RECH)
RCHZONE = ones(gr.Ny,gr.Nx); RCHZONE(:,gr.xm>gr.xGr(7))=2;

% parameters used here, for the two zones in RCHZONE
R=[par.prch1,par.prch2]*1e-8;

RECH=R(RCHZONE);  % efficient assignment of recharge values to zones

%% WELLS

wel=[  % well locations and extraction
     3	5	11  -5   % calibration parameter
     2	4	6   -5   % calibration parameter
     2	6	12  -5   % calibration parameter
     1	9	8   -5
     1	9	10  -5
     1	9	12  -5
     1	9	14  -5
     1	11	8   -5
     1	11	10  -5
     1	11	12  -5
     1	11	14  -5
     1	13	8   -5
     1	13	10  -5
     1	13	12  -5
     1	13	14  -5
];

%% DRAINS

% parameters used here as well
drn=[               %% drain locations
    1	8	2	0	1*par.pdrn  % conductance to be calibrated
    1	8	3	0	1*par.pdrn  % conductance to be calibrated
    1	8	4	10	1
    1	8	5	20	1
    1	8	6	30	1
    1	8	7	50	1
    1	8	8	70	1
    1	8	9	90	1
    1	8	10	100	1
    ];

iper=1;

%% Boundary conditions
WEL=[iper*ones(size(wel,1)) wel];
DRN=[iper*ones(size(drn,1)) drn];

save underneath IObs  % needed in mf_analyze to extract observations

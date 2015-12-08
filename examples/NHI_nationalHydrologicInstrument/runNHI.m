%% runNHI.m
% This file runs
% THis
% that

% NHI  % this runs mfLab/examples/NHI/NHI.m

%%
% which converts the stored ASCII data files in
% _mfLab/exmples/NHI/NHIascii_
% into mata files in
% _mfLab/examples/NHI/NHIdata_

%% NHI
% NHI is the "National Hydrological Instrument", the Dutch national
% groundwater-surface water model as published in
% _www.NHI.nu/bibliotheek_
%
% The exact locations (URLS filenames etc, can be found in the workbook
% NHI.xls in the current directory.
% See the readme's in the NHIascii and NHIzipfiles directory

%% Obtaining the NHI data set
% First thing to do after downloading the NHIzipfiles from the site
% _www.NHI.nu_
% is to run NHI.m, which reads the unpacked ASCII files from NHI
% from ./NHIzipfiles directory and converts them to
% mat files in the directory ./NHIdata.
%
% By doing this, the total volume of NHI files reduces from about 2GB
% to about 250 MB, while also reading data is about hundred times faster
% then reading the ASCII files directly or any known alternative.
% This way, we can keep the entire model in this data directory with
% while avoiding performance problems.

%% Obtaining submodels
% Make a new directory and name it properly. Copy the files
% _mf_adapt.m_ (sometimes _mf_build.m_)
% _mf_analyze.m_
% and the exel workbook
% to it and rename the workbook to reflect your particular model (perhaps
% giving it the same name as your model directory.
% In mf_adapt.m change the basenmae to this worksheet name. All model files
% will obtain this name, but with different extensions.
%
% To obtain a submodel of the NHI see one of the examples in this directory.
% These examples determine their center location from a Google
% Earth kml file containing a placmark (location). To get one start Google
% Earth, navigate to your area and put a placemark using a suitable name. It
% appears in the list next to the map. Select it, right-click, and save it
% as a kml file in the directory of your model.
% The mf_adapt.m file will extract the coordinates of the placemark and use
% them as the center of your model.
% in mf_adapt.m adapt the size of the model.
% Change anything else necessary in mf_adaptm, the Excel workbook and
% mf_analyze.m to obtain a correct submodel.

%% The examples
% The current directory contains a number of subdirectories each with an
% example named after its approximate location. Each one
% extracts a submodel from NHI and simulates it.
%
% To run such a model type mf_setup and when ready mf_analyze to see the
% results in 3D.

%% The examples
%% NL-West-East
% NL-West-East is a cross section through the center of the Netherlands
% from the North Sea to the German border, with the streamlines drawn in it.
% Only the aquitards are colored.

%% Amsterdam
% Amsterdam is a model to be used to anlayze thermal heat storage below Amsterdam.

%% Lexmond
% Lexmond is just a model in the center of the Netherlands of substantial
% size yet finishing in 60 seconds (steady state). It allows 3D viewing of
% a complex subsurface and also viewing of layer projections.
%
% The mf_analyze of this model also provides an extended number of figures
% of heads, transmissivities, and more.

%% Viewer
% There is much more that you can do with the 3D gui that its mf_analuze.m
% automatically invokes. You can launch the 3D viewer many times. It has
% to be deleted by hand, because the command _close_ _all_ has no
% effect on it.
%
% It is sometimes convenient to have more than one viewer on at the same
% time to show different variables simultaneously.

%%
% TO 120809
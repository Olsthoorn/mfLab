%% Creating Scenarios 
% This file will create scenarios for the model. First you can choose the
% number of runs for the model. It will stop automaticly if dS / dt < ... 
% [Mm3 / y]. After choosing the # of model runs you can select your
% basic scenario. Here it will be selected which features will be on or
% off. After the feature selecting you will need to give values for the
% additional information (default values are presented). You can select if
% you want more or less rainfall / recharge / pumping / infiltration. After
% the additional information, you can select if you would like to see plots
% of the model results. This can be done by false or true. Explained is
% provided to see what will be turned on. After all this information, it 
% will start running the model. All the information is read in from an
% excel sheet (JorfData.xls scenarios).
%
% This file needs to be used in combination with mf_setup, mf_adapt and 
% mf_analyze. Also for creating the waterbalance, Koen_waterbalance needs 
% to be used.  
%
% KG 171014

clear variables
close all
clc

%% del previous information
% these files will be created again in model run number 2 (if selected)
dos('del Waterbalance.mat');% here you delete the previous waterbalance
dos('del INITIALHD.mat');   % here you delete the previous storage (heads)
dos('del Scenarios.mat');   % here you delete the previous scenarios

%% Reading the scenarios
% Here you can choose your scenario. The scenario will be read for the file
% JorfData.xls. The file will be read, and it will find the values for the
% selected scenarios. The scenarios will be explained. The first 7 
% scenarios are with the default features. After these 7 basic scenarios 
% more detailed scenarios are given:

% Scenario 1:
    % This is the natural scenario. All the recharges are on, but there is
    % no human influence (no drains and pumps). The storage is made with
    % the default values, but those can be adjusted here in the scenarios.
    % The initial loaded storage is NaturalStorage.
% Scenario 2:
    % This is the human influence scenario. All features are turned on
    % (except Piezoms). The loaded storage will be the storage with the
    % default values (under Additional Information) and with all the
    % features on. The name of the loaded storage is HumanInfluenceStorage.    
% Scenario 3:
    % This is the natural scenario with floods. All the recharges are on, 
    % but there is no extraction (no drains and pumps). There are floodings 
    % in this scenario. The storage is made with the default values, but 
    % those can be adjusted here in the scenarios. The initial loaded 
    % storage is NaturalStorageFloods.
% Scenario 4:
    % This is the Khettara scenario (Full length Khettara). Everything is on 
    % except the pumping and piezoms. Additional information values are on
    % default. The storage is named KhettaraStorageFloods.
% Scenario 5:
    % This one is also Khettara storage. However, floods are turned off
    % here. Additional information values are on default. The storage is
    % named KhettaraStorage.
% Scenario 6:
    % This is the Pumping scenario. Everything is on except khettaras
    % (drains) and piezoms. Additional information values are on default.
    % The name of the storage is PumpingStorage.   
% Scenario 7:
    % This is the HumanInfluence scenario, but here the pumping in Oukhit
    % is included. Additional information is on default. the name of the
    % storage is HumanInfluenceStorageOukhit
    
% Scenario 8:
    % This scenario is based on HumanInfluenceStorage. However, the
    % recharge in this scenario will be based on the default precipitation
    % but is multiplied with a sinus function. This is to create random
    % precipitation to simulate dry and wet periods over the years.
% Scenario 9: 
    % Khettara history, how did the khettaras develop over time? The start is
    % with the natural storage (without floods) and khettaras are on.
    % However, the khettara length is limited. Over time the khettaras will
    % grow. The purpose is to see whether the khettaras became longer to
    % extract more water, or to see if they fell dry because of aquifer
    % depletion. Initial parameter values are on default. S.KhetLength will
    % determine the growth of the khettaras for the timestep.
% Scenario 10:
    % Random numbers: here the recharge (rain,floods,riverinfiltration) is
    % changed. This is to model dry or wet periods. It is based on a sinus.
    % The return period can be influenced so different results (longer or
    % shorter dry or wet periods) are different. In this scenario Pumping
    % and floods are turned off, and khettara length is turned on. This is
    % to see what happened in the past with dry and wet periods. The
    % starting storage is NaturalStorage. Note that this scenario is
    % exactly the same as scenario 8, but with random numbers turned on.
% Scenario 11: 
    % Pumping effects, how does pumping affect khettara discharge over
    % time? Starting with Khettara storage Floods, pumping will be turned 
    % on to see how this effects the discharge of khettaras. Initial 
    % parameter values are on default. Floods are turned on in this 
    % scenario. 
% Scenario 12:
    % Oukhit will start pumping (more) water from the aquifer. How will
    % this influence the area? Values are on default and the starting
    % storage is HumanInfluenceStorage. TimePump = true but only for the
    % Oukhit area. The pumping rate in Fezna Jorf is on a default (or non
    % changeble) rate.
% Scenario 13:
    % This scenario everything is turned on, except khettara length and
    % pumping time and piezoms. This is to see what happens to the system
    % with different precipitation rates (wet and dry periods). the
    % starting storage is HumanInfluenceStorageOukhit.
% Scenario 14:
    % In this scenario the super khettaras will be turned on. But only the 
    % ones in Jorf and Hannabou. They will replace the khettaras which are
    % currently there. Starting storage is HumanInfluenceStorage and values
    % are on default. Everything is turned on.
% Scenario 15:
    % In this scenario the super khettara of Jorf is turned on. The old
    % khettaras remain functional and also the other features are turned
    % on. values are on default and the starting storage is
    % HumanInfluenceStorage.
% Scenario 16:
    % In this scenario additional infiltration upstream of Fezna is turned
    % on. Storage is HumanInfluenceStorage and values are on default. This
    % is to see what happens when additional water is infiltration from the
    % river into the aquifer.

% Scenario 17:
    % This scenario is empty. Only the head boundaries are turned on. In
    % the excel sheet JorfData.xls in the map XLSData values in this
    % scenario can be changed, so you can create your own scenario. Values
    % of the features can be changed in this matlab file.
    
% Choose the scenario:
Scenario = 15;

XLSData = ['.' filesep 'XLSData'  filesep]; % EXCEL files directory
MATFiles= ['.' filesep 'MATFiles' filesep]; % Storage files directory
MFiles  = ['.' filesep 'MFiles'   filesep]; % Storage files directory

addpath(MATFiles); % Add to path so they can be located
addpath(MFiles);   % Add to path so they can be located

% Read in the Scenarios sheet in the excel file JorfData.xls
[ScenarioVAL,ScenarioNAM] = xlsread([XLSData,'JorfData.xls'],'Scenarios');

% Get the scenario values (true of false) and the values for the parameters
[S] = GetScenarioVals(Scenario,ScenarioNAM,ScenarioVAL);

%% Save scenario to load it into the model.
save Scenarios S

%% Run the model
Run_Model(S);

fprintf('\n\n.. Finished\n\n\n\n');
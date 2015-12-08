function [S] = GetScenarioVals(Scenario,ScenarioNAM,ScenarioVAL)
% This function will get the values for the specified scenario. It will
% also get information for the storage and it will try to find a previous
% waterbalance associated with the storage.
%
% KG 141103

%% Get the row number of the features.
    NumberOfRuns    = strmatchi('NumberOfRuns',ScenarioNAM);
    MaxStorCh       = strmatchi('MaxStorCh',ScenarioNAM);

    Pump            = strmatchi('Pumping',ScenarioNAM);
    Floods          = strmatchi('Floods',ScenarioNAM);
    Khettaras       = strmatchi('Khettaras',ScenarioNAM);
    River           = strmatchi('River',ScenarioNAM);
    Recharge        = strmatchi('Recharge',ScenarioNAM);
    EVT             = strmatchi('EVT',ScenarioNAM);
    ConstantH       = strmatchi('ConstantH',ScenarioNAM);
    GradientBS      = strmatchi('GradientBS',ScenarioNAM);
    Piezom          = strmatchi('Piezom',ScenarioNAM);
    UseLH           = strmatchi('UseLH',ScenarioNAM);
    Storage         = strmatchi('Storage',ScenarioNAM);
    KhetLength      = strmatchi('KhetLength',ScenarioNAM);
    TimePump        = strmatchi('TimePump',ScenarioNAM);
    Random          = strmatchi('Random',ScenarioNAM);
    Oukhit          = strmatchi('Oukhit',ScenarioNAM);
    SuperKhetFJ     = strmatchi('SuperKhetFJ',ScenarioNAM);
    SuperKhetHa     = strmatchi('SuperKhetHa',ScenarioNAM);
    SuperKhetFezna  = strmatchi('SuperKhetFezna',ScenarioNAM);
    HoldKhettaras   = strmatchi('HoldKhettaras',ScenarioNAM);
    SandInf         = strmatchi('SandInf',ScenarioNAM);
    
    % Additional parameters
    LayerInf        = strmatchi('Layerinf',ScenarioNAM);
    Dsurf           = strmatchi('Dsurf',ScenarioNAM);
    ExDepth         = strmatchi('ExDepth',ScenarioNAM);
    InfFact         = strmatchi('InfFact',ScenarioNAM);
    InitialLoss     = strmatchi('InitialLoss',ScenarioNAM);
    RunOff          = strmatchi('RunOff',ScenarioNAM);
    sy              = strmatchi('sy',ScenarioNAM);
    
    PumpVol         = strmatchi('PumpVolFJ',ScenarioNAM);
    PumpVolCh       = strmatchi('PumpVolCh',ScenarioNAM);
    PumpVolMax      = strmatchi('PumpVolMax',ScenarioNAM);
    YearPump        = strmatchi('YearPump',ScenarioNAM);
    StartYPump      = strmatchi('StartYPump',ScenarioNAM);
    PumpVolOuk      = strmatchi('PumpVolOuk',ScenarioNAM);
    PumpOukhitMax   = strmatchi('PumpOukhitMax',ScenarioNAM);
    PumpVolChOuk    = strmatchi('PumpVolCOuk',ScenarioNAM);
    
    MaxLengthKhet   = strmatchi('MaxLengthKhet',ScenarioNAM);
    Length          = strmatchi('Length',ScenarioNAM);
    YearKhet        = strmatchi('YearKhet',ScenarioNAM);
    StartYKhet      = strmatchi('StartYKhet',ScenarioNAM);
    EindYKhet       = strmatchi('EindYKhet',ScenarioNAM);
    NumberOfKhetRuns = strmatchi('NumberOfKhetRuns',ScenarioNAM);
    
    PercentFlooding = strmatchi('PercentFlooding',ScenarioNAM);
    PercentRecharge = strmatchi('PercentRecharge',ScenarioNAM);
    PercentRiver    = strmatchi('PercentRiver',ScenarioNAM);
    PercentEVT      = strmatchi('PercentEVT',ScenarioNAM);
    
    FloodingFJ   = strmatchi('FloodingFJ',ScenarioNAM);
    FloodingHa   = strmatchi('FloodingHa',ScenarioNAM);
    FloodingSand = strmatchi('FloodingSand',ScenarioNAM);   

    % Additional plots
    LayerThick      = strmatchi('LayerThick',ScenarioNAM);
    Control         = strmatchi('Control',ScenarioNAM);
    ContourDEM      = strmatchi('ContourDEM',ScenarioNAM);
    LayerDEM        = strmatchi('LayerDEM',ScenarioNAM);
    ProfilePump     = strmatchi('ProfilePump',ScenarioNAM);
    Dview           = strmatchi('Dview',ScenarioNAM);
    TotalTrans      = strmatchi('TotalTrans',ScenarioNAM);
    SpyCatch        = strmatchi('SpyCatch',ScenarioNAM);
    Head            = strmatchi('Head',ScenarioNAM);
    Map             = strmatchi('Map',ScenarioNAM);
    KhetPlot        = strmatchi('KhetPlot',ScenarioNAM);
    
    Cross           = strmatchi('CrossAll',ScenarioNAM);
    CrossKhet       = strmatchi('CrossKhet',ScenarioNAM);
    Qkhettaras      = strmatchi('Qkhettaras',ScenarioNAM);
    KhetIntersect   = strmatchi('KhetIntersect',ScenarioNAM);

%% get starting info model:

S.NumberOfRuns = ScenarioVAL(NumberOfRuns,Scenario);
S.MaxStorCh    = ScenarioVAL(MaxStorCh,Scenario);

%% Turning features on or off:
% Select which features you want to run in the model. true (1) means it 
% will be on, false (0) means it will be off.

% Will turn pumping (-pumping infiltration) on or off
    S.Pumping      = ScenarioVAL(Pump,Scenario); 
% Will turn flooding on or off
    S.Floods       = ScenarioVAL(Floods,Scenario); 
% Will turn khettara drainage / infiltration on or off
    S.Khettaras    = ScenarioVAL(Khettaras,Scenario); 
% Will turn on or off river infiltration
    S.River        = ScenarioVAL(River,Scenario); 
% Will turn on / off recharge in the area
    S.Recharge     = ScenarioVAL(Recharge,Scenario); 
% Will turn on / off Evapotransperation
    S.EVT          = ScenarioVAL(EVT,Scenario); 
% Turns on / off Constant Head boundaries
    S.ConstantH    = ScenarioVAL(ConstantH,Scenario); 
% Turns on / off Gradient boundary at south boundary
    S.GradientBS   = ScenarioVAL(GradientBS,Scenario);
% Will turn on / off the piezometers and plots of piezometers
    S.Piezom       = ScenarioVAL(Piezom,Scenario); 
% Use last heads generated to be used again (has large effect on storage)
    S.UseLH        = ScenarioVAL(UseLH,Scenario); 
% Make the khettara length changeable
    S.KhetLength   = ScenarioVAL(KhetLength,Scenario);
% Makes Pumping variable over time
    S.TimePump     = ScenarioVAL(TimePump,Scenario);
% Will use random rainfall (uses a sin function)
    S.Random       = ScenarioVAL(Random,Scenario);
% Will turn on pumping around Oukhit, Volume needs to be defined.
    S.Oukhit       = ScenarioVAL(Oukhit,Scenario);
% Will turn on super khettara in Jorf    
    S.SuperKhetFJ  = ScenarioVAL(SuperKhetFJ,Scenario);
% Will turn on super khettara in Fezna Jorf
    S.SuperKhetFezna= ScenarioVAL(SuperKhetFezna,Scenario);
% will turn on super khettara in Hannabou
    S.SuperKhetHa  = ScenarioVAL(SuperKhetHa,Scenario);
% Will keep other khettaras with the super khettara in Jorf
    S.HoldKhettaras= ScenarioVAL(HoldKhettaras,Scenario);
% will turn on sand infiltration parts in upstream area
    S.SandInf      = ScenarioVAL(SandInf,Scenario);

%% Selecting values Additional information (parameters)

S.LayerInf      = ScenarioVAL(LayerInf,Scenario);
S.DSurf         = ScenarioVAL(Dsurf,Scenario);
S.ExDepth       = ScenarioVAL(ExDepth,Scenario);
S.InfFact       = ScenarioVAL(InfFact,Scenario);
S.InitialLoss   = ScenarioVAL(InitialLoss,Scenario);
S.RunOff        = ScenarioVAL(RunOff,Scenario);
S.sy            = ScenarioVAL(sy,Scenario);

S.PumpVol       = ScenarioVAL(PumpVol,Scenario);
S.PumpVolCh     = ScenarioVAL(PumpVolCh,Scenario);
S.PumpVolMax    = ScenarioVAL(PumpVolMax,Scenario);
S.YearPump      = ScenarioVAL(YearPump,Scenario);
S.StartYPump    = ScenarioVAL(StartYPump,Scenario);
S.PumpVolOuk    = ScenarioVAL(PumpVolOuk,Scenario);
S.PumpOukhitMax = ScenarioVAL(PumpOukhitMax,Scenario);
S.PumpVolChOuk  = ScenarioVAL(PumpVolChOuk,Scenario);

S.MaxLengthKhet = ScenarioVAL(MaxLengthKhet,Scenario);
S.Length        = ScenarioVAL(Length,Scenario);
S.YearKhet      = ScenarioVAL(YearKhet,Scenario);
S.StartYKhet    = ScenarioVAL(StartYKhet,Scenario);
S.EindYKhet     = ScenarioVAL(EindYKhet,Scenario);
S.NumberOfKhetRuns = ScenarioVAL(NumberOfKhetRuns,Scenario);

S.PercentFlooding = ScenarioVAL(PercentFlooding,Scenario);
S.PercentRecharge = ScenarioVAL(PercentRecharge,Scenario);
S.PercentRiver    = ScenarioVAL(PercentRiver,Scenario);
S.PercentEVT      = ScenarioVAL(PercentEVT,Scenario);

S.FloodingFJ      = ScenarioVAL(FloodingFJ,Scenario);
S.FloodingHa      = ScenarioVAL(FloodingHa,Scenario);
S.FloodingSand    = ScenarioVAL(FloodingSand,Scenario);

S.LayerThick      = ScenarioVAL(LayerThick,Scenario);
S.Control         = ScenarioVAL(Control,Scenario);
S.ContourDEM      = ScenarioVAL(ContourDEM,Scenario);
S.LayerDEM        = ScenarioVAL(LayerDEM,Scenario);
S.ProfilePump     = ScenarioVAL(ProfilePump,Scenario);
S.DView           = ScenarioVAL(Dview,Scenario);
S.TotalTrans      = ScenarioVAL(TotalTrans,Scenario);
S.SpyCatch        = ScenarioVAL(SpyCatch,Scenario);
S.Head            = ScenarioVAL(Head,Scenario);
S.Map             = ScenarioVAL(Map,Scenario);
S.KhetPlot        = ScenarioVAL(KhetPlot,Scenario);

S.Cross           = ScenarioVAL(Cross,Scenario);
S.CrossKhet       = ScenarioVAL(CrossKhet,Scenario);
S.QKhettaras      = ScenarioVAL(Qkhettaras,Scenario);
S.KhetIntersect   = ScenarioVAL(KhetIntersect,Scenario);

%% Get the rainfall data from 2000 - 2011
% Use this rainfall - initial loss to determine the factor for recharge
% with flooding scenarios.

basename= 'Jorf';
[PERnams,PERvals,NPER] = getPeriods(basename,'PER',{'PERLEN','Precip'});
precip = PERvals(:,strmatchi('precip',PERnams));
perlen = PERvals(:,strmatchi('PERLEN',PERnams));
S.RainRech = max(precip,0) / NPER / perlen(1);

%% Initial Storage
% The model will start with starting heads and will the continue from here
% on. The starting heads will define the initial storage, therefore the
% initial storage should be changeable. In the excel sheet with the
% scenarios the number for the storage is defined. An initial storage can
% be usefull to reduce computational time to get to an equilibrium. Also,
% the waterbalance is added to the storage. This will further reduce
% computational time due to feedback features.

% load a previous storage, if none is selected the model will create one.
n = ScenarioVAL(Storage,Scenario); % here you can select the scenario, look below to see the storages

switch n
    case 1
        % Natural storage (no floods)
        try
            load NaturalStorage
            fprintf('... Loading Natural Storage');
        catch
            fprintf('.. Natural Storage not found');
        end
    case 2
        try
            load HumanInfluenceStorage                 
            fprintf('... Loading Human Influence Storage\n');
        catch
            fprintf('.. HumanInfluenceStorage not found');
        end
    case 3
        try
            load NaturalStorageFloods                 
            fprintf('... Loading Natural Storage with floods\n');
        catch
            fprintf('.. NaturalStorageFloods  not found');
        end
    case 4
        try
            load KhettaraStorageFloods
            fprintf('... Loading Khettara Storage with floods\n');
        catch
            fprintf('.. Khettara Storage with Floods not found');
        end
    case 5
        try
            load KhettaraStorage
            fprintf('... Loading Khettara Storage\n');
        catch
            fprintf('.. Khettara Storage not found');
        end
    case 6
        try
            load PumpingStorage         
            fprintf('... Loading Pumping Storage\n'); 
         catch
            fprintf('.. Pumping Storage not found');
        end
    case 7
        try
            load HumanInfluenceStorageOukhit
            fprintf('... Loading Human Influence Storage Oukhit\n'); 
        catch
            fprintf('.. HumanInfluence Storage Oukhit not found');
        end
end

% Save it to InitialHD, if this doesn't work the model will create Initial
% heads to start with. (if INITIALHD is not deleted earlier on).
try
    save INITIALHD STRTHD xc yc
catch 
    fprintf('No initial storage (heads) found\n');
end

% Also try to save waterbalance (If available in the storage). If not it
% will be created in the model.
try
    save Waterbalance W
catch
    fprintf('No previous waterbalance found\n');
end

%% Random precipitation
if S.Random
    % will create 'random' precipitation. it will be in a % of the data
    % precipitation. It will be plotted here.
    if S.KhetLength
        n11 = (S.EindYKhet - S.StartYKhet)/11;
        S.StartYKhet = S.EindYKhet - ceil(n11) * 11;
        t =  (S.StartYKhet:11:S.EindYKhet)';  % make sure EindYKhet is included
    else
        t = (1:S.NumberOfRuns)';
    end

    if false
        T = [19 20 21]'; % Return periods

        % Select the period. Let it be dependend on Khettara building or the
        % number of runs.

        % Calculate the factor y.
        % additional factor to get more extreme values 
        S.y = sum(sin(2 * pi .* (1./T)*t + 2 * pi) / (numel(T))) + 1;
    else
        load('precipFac');
        I    = precipFac(:,1)>S.StartYKhet-11; %#ok
        yFac = precipFac(I,end);
        S.y  = reshape(yFac,[11,numel(yFac)/11]);
        % 11xperiods ieach column is an eleven year period, it has a
        % differnt factor for each year in the given period, exactly
        % according to Till and Guiot (1990)
    end
    
    figure;
    plot(t,S.y); ylim([0 2]);
end
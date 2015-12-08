%% Added by Koen: Waterbalance.
% this script will calculated the waterbalance in the area. It will also
% create waterbalances for local areas such as Fezna & Jorf and Hannabou.
% It will create a water balance for all the features. If features are
% turned off in scenarios the water balance will show an outflow (or
% inflow) of 0.
%
% KG 171014

load Scenarios     % load the scenarios
load name          % get name of model
load(basename)     % load model data save in mf_setup
load underneath    % get DEM etc, saved at end of mf_adapt

Ayear = 365.24;    % 1 year is 365.24 days

%% Save the last head to be used as initial heads in a subsequent run
% Here the Last Head will be saved. It will provide the starting head for
% the next run.
H = readDat([basename '.HDS']); % get the heads
B = readBud([basename '.BGT'],'t',[H.totim]); % read the budget file

if S.UseLH
    STRTHD   = H(end).values(:,:,1);   % #ok
    strthd   = gr.Z(:,:,end) + 0.5 * abs(diff(gr.Z(:,:,([1 end])),1,3));
    STRTHD(isnan(STRTHD)) = strthd(isnan(STRTHD));
    xc = gr.xc;
    yc = gr.yc;
    save INITIALHD STRTHD xc yc
end   
        
%% In and out groundwater fluxes North, South and anti atlas
% Here the balance for the North and South boundary will be calculated and
% presented. First you locate the places of the boundaries (by looking in
% headboundaries. here the Idx is present). Then you will find the location
% of these cells in the budget file. The values will be selected. The sum
% of the values over all the timesteps will be the balance for the
% boundaries.

if ~S.GradientBS
    Bnorth = strmatchi('north',{headBoundaries.name}); % find location of the name
    Bsouth = strmatchi('south',{headBoundaries.name}); % find location of the name
    TotN = zeros(size(B)); TotS = zeros(size(B));      % create zero matrix

    for iB=1:length(B)
        Head = strmatchi('constanthead',B(iB).label); % Find the correct label
        Head = B(iB).term{Head};                      % find the label values

        HeadN = Head(headBoundaries(Bnorth).Idx); % only for selected cells of boundary
        HeadS = Head(headBoundaries(Bsouth).Idx);

        TotN(iB) = sum(HeadN); % sum of the cells
        TotS(iB) = sum(HeadS);
    end

    W.QNorth    = TotN * Ayear / 1000000; % [Mm3 / y] at every time step
    W.IN.QNorth = mean(W.QNorth);         % average [Mm3 / y]
    W.QSouth 	= TotS * Ayear / 1000000; % [Mm3 / y] at every time step
    W.OUT.QSouth= mean(W.QSouth);         % average [Mm3 / y]

        fprintf('Total inflow at North Boundary:\n');
        fprintf('Qtot North [Mm3/y]-->'); 
            fprintf('%12.2f',W.IN.QNorth); fprintf('\n\n');

        fprintf('Total outflow at South Boundary:\n');
        fprintf('Qtot South [Mm3/y]-->'); 
            fprintf('%12.2f',W.OUT.QSouth); fprintf('\n\n');
else
    Bnorth = strmatchi('north',{headBoundaries.name}); % find location of the name
    TotN = zeros(size(B));
    
    for iB=1:length(B)
        Head = strmatchi('constanthead',B(iB).label); % Find the correct label
        Head = B(iB).term{Head};                      % find the label values
        HeadN = Head(headBoundaries(Bnorth).Idx);
        TotN(iB) = sum(HeadN); % sum of the cells
    end
    
    W.QNorth    = TotN * Ayear / 1000000; % [Mm3 / y] at every time step
    W.IN.QNorth = mean(W.QNorth);         % average [Mm3 / y]
    
        fprintf('Total inflow at North Boundary:\n');
        fprintf('Qtot North [Mm3/y]-->'); 
            fprintf('%12.2f',W.IN.QNorth); fprintf('\n\n');
end

%% Recharge
% This part will calculate how much recharge there will be in the model.
% First it will look at the correct label in the budget file, then it will
% get the corresponding values from the budget file. The recharge for every
% time step will be calculated.

Recharge = zeros(size(B)); % Create the zero matrix

for iR=1:length(B)
    Nrech    = strmatchi('recharge',B(iR).label);
    rech        = B(iR).term{Nrech};
    Recharge(iR)= sum(rech(:)); % Sum at every time step
end

W.QRecharge = Recharge * Ayear / 1000000; % [Mm3 / y] at every time step
W.IN.QRecharge = mean(W.QRecharge);       % Average [Mm3 / y]

    fprintf('Total recharge in the Area:\n');
    fprintf('Recharge [Mm3/y] -->'); 
        fprintf('%12.2f',W.IN.QRecharge); fprintf('\n\n');

%% River infiltration
% The river has been divided into two parts. First you will get the
% location of the rivers from the riversQ file. Then you will look up the
% values in the budget file. Not every time step will have river
% infiltration. This is dependend on precipitation. If there is no
% infiltration the time step value will be 0 (for both rivers)
River1 = riversQ(1).Idx;
River2 = riversQ(2).Idx;

for iR = length(B):-1:1
    RivL = strmatchi('riverleakage',B(iR).label);
    River = B(iR).term{RivL};

    QRiv1(iR) = sum(River(River1));
    QRiv2(iR) = sum(River(River2));
end

W.QRiv1 = mean(QRiv1) * Ayear / 1000000;    % [Mm3 / y] at every time step
W.QRiv2 = mean(QRiv2) * Ayear / 1000000;    % [Mm3 / y] at every time step
W.IN.QRiver = mean(W.QRiv1) + mean(W.QRiv2);% average [Mm3 / y]

    fprintf('Total infiltration from river in the Area:\n');
    fprintf('River infiltration [Mm3/y] -->'); 
        fprintf('%12.2f',W.IN.QRiver); fprintf('\n\n');

%% Storage
% Here the storage from the budget file will be read. It is devided in
% negative and positive storage. If the storage is positive; This means
% there is water flowing out of storage. If the storage is negative water
% will be flowing into storage. The sum of the storage should be positive -
% negative. The sum of the storages is the change in storage.

StorageP = zeros(size(B)); StorageN = zeros(size(B)); som = zeros(size(B));

for iS = 1:length(B)
    St = strmatchi('storage',B(iS).label);
    Stor = B(iS).term{St};
    som(iS) = sum(Stor(:));
    
    Negative = Stor < 0; % select negative values
    Positive = Stor > 0; % select positive values
    
    StorageP(iS) = sum(Stor(Positive));
    StorageN(iS) = sum(Stor(Negative));
end

W.StorageP     = StorageP * Ayear / 1000000; % [Mm3 / y] at every time step
W.StorageN     = StorageN * Ayear / 1000000; % [Mm3 / y] at every time step
W.StorageSum   = som * Ayear / 1000000;      % [Mm3 / y] at every time step
W.IN.StorageP  = mean(W.StorageP);           % average [Mm3 / y]
W.OUT.StorageN = mean(W.StorageN);           % average [Mm3 / y]

    fprintf('Change in storage:\n');
    fprintf('Storage Positive [Mm3/y] -->'); 
        fprintf('%12.2f',W.IN.StorageP); fprintf('\n');
    fprintf('Storage Negative [Mm3/y] -->'); 
        fprintf('%12.2f',W.OUT.StorageN); fprintf('\n');
    fprintf('Change Storage   [Mm3/y] -->'); 
        fprintf('%12.2f',mean(W.StorageSum)); fprintf('\n\n');

        
%% Khettaras
% Here the extraction of the khettaras will be calculated. Also, the
% khettara extraction for Fezna & Jorf and for Hannabou will calculated
% seperately. This will provide feedback for the infiltration for the next
% model run. Because some khettaras drain from the same cell, the unique
% location cells must be found first. This is done in the loop. The name of
% the area where the khettara drain is found in the excel sheet
% JorfData.xls. The withdrawel from Fezna and Jorf + Hannabou should equal
% to the sum of the khettaras.

QKhettara = zeros(size(B));   QKhettaraFJ = zeros(size(B)); 
QKhettaraHa = zeros(size(B)); QKhettaraF  = zeros(size(B)); 

% for seeking the unique location values, it will also look at cell number
% 1 for khettara drainage, but here the ibound = 0 and it is far away for
% where the khettaras actually are. Therefore, you can use 1 for the
% starting value.
HaIdx = 1; FJIdx = 1; IdxKhet = 1; FIdx = 1;

% It could be that S.GradientBS is selected. If so, the south boundary has
% become a drain. Therefor, the location needs to save in order to find the
% correct idx. This part is always true, but in the loop it will only be
% done if S.GradientBS = true.
SouthB = strmatchi('s',{headBoundaries.name});
QSouthB = zeros([length(B) length(SouthB)]);

for iB=1:length(B)
    drains = strmatchi('drains',B(iB).label); % find location of drains in header
    Khet = B(iB).term{drains};                % get the balance
    
    
    if S.GradientBS
        for iLay = 1:length(SouthB)
            QSouthB(iB,iLay) = sum(Khet(headBoundaries(SouthB(iLay)).Idx));
        end
        QKhettara(iB) = sum(Khet(:)) - sum(QSouthB(iB,:));
    else
        QKhettara(iB) = sum(Khet(:));  % make a sum of the balance
    end
    
    FindIdxKhet = find(Khet);
    IdxKhet = unique([IdxKhet FindIdxKhet']); % used for making plot of 
                                              % draining cells

    for iK=1:length(khettaras)
        % This loop is to seek in which infiltration area the khettaras
        % drain.
        IArea = strmatchi('Infiltratiegebied naam',khettaras(iK).vertexTxtHdr);
        NArea = khettaras(iK).vertexTxt(1,IArea);

        if strcmp('Hannabou',NArea)
            % This will get all the cell locations from the khettaras which
            % drain towards Hannabou, the same is done for the khettaras
            % which drain towards Fezna and Jorf
            HaIdx = [HaIdx khettaras(iK).Idx]; 
        elseif strcmp('Fezna Jorf',NArea)
            FJIdx = [FJIdx khettaras(iK).Idx];
        elseif strcmp('fezna',NArea)
            FIdx = [FIdx khettaras(iK).Idx];
        else
            fprintf('\n  ERROR!!!!!\n');
        end
    end
    
    % Here the balance for the unique values are calculated for every time
    % step in the model.
    QKhettaraFJ(iB) = sum(Khet(unique(FJIdx)));
    QKhettaraHa(iB) = sum(Khet(unique(HaIdx)));
    QKhettaraF(iB)  = sum(Khet(unique(FIdx)));

    % This is to clear the values it contains for the next loop.
    HaIdx = 1; FJIdx = 1; FIdx = 1;
end
    
W.QKhettara     = QKhettara * Ayear / 1000000;   % [Mm3 / y] at every time step
W.QKhettaraFJ   = QKhettaraFJ * Ayear / 1000000; % [Mm3 / y] at every time step
W.QKhettaraHa   = QKhettaraHa * Ayear / 1000000; % [Mm3 / y] at every time step
W.QKhettaraF    = QKhettaraF * Ayear / 1000000;  % [Mm3 / y] at every time step
W.OUT.QKhettara = mean(W.QKhettara);             % [Mm3 / y]

    fprintf('Total withdrawal by khettaras:\n');
    fprintf('Qtot khettaras   [Mm3/y]-->'); 
        fprintf('%12.2f',W.OUT.QKhettara); fprintf('\n');
    fprintf('Q Fezna Jorf     [Mm3/y]-->'); 
        fprintf('%12.2f',mean(W.QKhettaraFJ)); fprintf('\n');
    fprintf('Q Fezna Hannabou [Mm3/y]-->'); 
        fprintf('%12.2f',mean(W.QKhettaraHa)); fprintf('\n'); 
    fprintf('Q Fezna          [Mm3/y]-->'); 
        fprintf('%12.2f',mean(W.QKhettaraF)); fprintf('\n'); 

if S.GradientBS
    W.QSouth        = mean(SouthB(:)) * Ayear / 1000000;
    W.OUT.QSouth    = W.QSouth;
    
    fprintf('Q South Boundary [Mm3/y]-->'); 
        fprintf('%12.2f',mean(W.QSouth)); fprintf('\n\n');         
end    
%% Pumping and infiltration
% Pumping and infiltration happens at the same place in the budget file. It
% happens in WELLS in the budget file. In WELLS infiltration and extraction
% can happen in the same cell. It is difficult to get the correct values
% from the files. If it happens in the same layer, a balance is found for
% the area, but constant rates (Pump infiltration and Flood infiltration)
% are substracted from this balance. The constant rate in the first loop is
% based on the initial values. In the next loop it is based on the
% waterbalance. This is done in mf_adapt.
% The model boundaries have been change in the process of making the model.
% The south boundary change, and therefor you need to compare the IBOUND
% with the Idx of Hannabou and only select the Idx which is in the active
% IBOUND.
% An infiltration factor is used to determine how much water will return to
% the model. about 60 - 70% of the water drained or pump will evaporate,
% the rest of the water will infiltrate. However, if the groundwatertable
% is close to the surface, more water will evaporate.
% The results will be presented later on, after the ET is calculated. This
% is to give a full balance of the local areas.
TotWell = zeros(size(B));

% Search for names related to pumping
    FJP     = strmatchi('FeznaJorfPump',{pumpAreas.name});   
    FJPI    = strmatchi('FeznaJorfPum-Inf',{pumpAreas.name}); 
    OHP     = strmatchi('Oukhit',{pumpAreas.name});   
    OHPI    = strmatchi('Oukhi-Inf',{pumpAreas.name});     
    BFJP    = zeros(size(B)); BFJPI = zeros(size(B));
    BOHP    = zeros(size(B)); BOHPI = zeros(size(B));
    IFluxP  = strmatchi('flux',pumpAreas.vertexHdr);

% Search for names related to Flood infiltration
[PERnams,PERvals,NPER] = getPeriods(basename,'PER',{'preCip'});
        
    HFI     = strmatchi('Hannabou-flood',{FloodAreas.name});      
    FJFI    = strmatchi('FeznaJorfflood-Inf',{FloodAreas.name}); 
    Sand    = strmatchi('SandInf',{FloodAreas.name});
    BHFI    = zeros(size(B)); BFJFI = zeros(size(B));
    BSI     = zeros([length(B) length(Sand)]);
    IFluxF  = strmatchi('flux',FloodAreas.vertexHdr);
    
if S.Floods
    FJFloodFactor = (S.FloodingFJ / sum(FloodAreas(FJFI).A) / 365.24 ) /...
        sum(S.RainRech);
    HaFloodFactor = (S.FloodingHa / sum(FloodAreas(HFI).A) / 365.24 ) /...
        sum(S.RainRech);
    if S.SandInf
        SandFloodFactor = (S.FloodingSand / (sum(FloodAreas(Sand(1)).A) + ...
            sum(FloodAreas(Sand(2)).A) + sum(FloodAreas(Sand(3)).A)) ...
            /Ayear )/sum(S.RainRech);
    else
        SandFloodFactor = 0;
    end
else
    FJFloodFactor = 0;
    HaFloodFactor = 0;
    SandFloodFactor = 0;
end

% Search for names related to khettara infiltration
    FJKI    = strmatchi('FeznaJorfKhet-Inf',{pumpAreas.name});  
    HKI     = strmatchi('HannabouKhet-Inf',{pumpAreas.name});
    FKI     = strmatchi('FeznaKhet-Inf',{pumpAreas.name}); 
    BHKI    = zeros(size(B)); BFJKI = zeros(size(B)); BFKI = zeros(size(B));
    IFluxK  = strmatchi('flux',pumpAreas.vertexHdr);

if S.LayerInf == 3.5 % if Wells (and infiltration + pump) is in same layer
    for iW = 1:length(B)
        wel = strmatchi('wells',B(iW).label);
        Well = B(iW).term{wel};

        % Pumping and pumping infiltration
        % Pumping rate = sum in the area - pump infiltration - flood
        % infiltration;
        BFJP(iW)  = sum(Well(pumpAreas(FJP).Idx) - mean(pumpAreas(FJP).A)*...
            pumpAreas(FJPI).vertex(1,IFluxP) - mean(FloodAreas(FJFI).A)*...
            FloodAreas(FJFI).vertex(1,IFluxF) - W.QKhettaraJ(iW));
        BFJPI(iW) = sum(pumpAreas(FJPI).A*pumpAreas(FJPI).vertex(1,IFluxP));
        
%         BOHP(iW) = sum(Well(pumpAreas(OHP).Idx)); ..
%         BOHPI(iW)= sum(Well(pumpAreas(OHPI).Idx));

        % Flood and flood infiltration
        BFJFI(iW) = PERvals(iW) * FJFloodFactor * S.InfFact * ...
            sum(FloodAreas(FJFI).A);
        BHFI(iW) = PERvals(iW) * HaFloodFactor * S.InfFact * ...
            numel(intersect(find(IBOUND~=0),FloodAreas(HFI).Idx)) *...
            mean(FloodAreas(HFI).A);
        
        for iS = 1:length(Sand)
            BSI(iW,iS) = Well(FloodAreas(Sand(iS)).Idx);
        end

        % Khettera and Khettera infiltration
        BFJKI(iW) = sum(Well(pumpAreas(FJKI).Idx));
        BHKI(iW)  = sum(Well(intersect(find(IBOUND~=0) ,...
            pumpAreas(HKI).Idx))) - BHFI(iW);
    end
else % if wells (the infiltration + pump) is located in different layers
    for iW = 1:length(B)
        wel = strmatchi('wells',B(iW).label);
        Well = B(iW).term{wel};
        
        TotWell(iW) = sum(B(iW).term{wel}(:));
        
        % Flood and flood infiltration
        BFJFI(iW) = PERvals(iW) * S.PercentFlooding / 30.5 * ...
            FJFloodFactor * S.InfFact * sum(FloodAreas(FJFI).A);
        BHFI(iW) = PERvals(iW) * S.PercentFlooding / 30.5 * ...
            HaFloodFactor * S.InfFact * numel(intersect(find(IBOUND~=0),...
            FloodAreas(HFI).Idx)) * mean(FloodAreas(HFI).A);
        
        for iS = 1:length(Sand)
            BSI(iW,iS) =  PERvals(iW) * S.PercentFlooding / 30.5 * ...
                SandFloodFactor * S.InfFact * sum(FloodAreas(Sand(iS)).A);
        end

        % Pumping and pumping infiltration
        BFJP(iW) = sum(Well(pumpAreas(FJP).Idx));
        BFJPI(iW)= sum(sum(pumpAreas(FJPI).A)*pumpAreas(FJPI).vertex(1,IFluxP));
        BFKI(iW) = sum(Well(pumpAreas(FKI).Idx)) - BFJPI(iW) - BFJFI(iW);
            
        BOHP(iW) = sum(Well(pumpAreas(OHP).Idx));
        BOHPI(iW)= sum(Well(pumpAreas(OHPI).Idx));
        
        % Khettera and Khettera infiltration
        BFJKI(iW) = sum(Well(pumpAreas(FJKI).Idx));
        BHKI(iW)  = sum(Well(intersect(find(IBOUND~=0) , ...
            pumpAreas(HKI).Idx))) - BHFI(iW);
    end 
end

% all in [Mm3 / y] per time step:
W.TotWell = sum(TotWell) / numel(B) * Ayear / 1000000; % total balance for control
W.QFJP  = BFJP * Ayear / 1000000;  % Fezna Jorf pump area balance (pump)
W.QFJPI = BFJPI * Ayear / 1000000; % Fezna Jorf pump area balance (pump inf)
W.QFKI  = BFKI * Ayear / 1000000;  % Jorf khettara inf (pump area)
W.QOHP  = BOHP * Ayear / 1000000;  % Oukhit pump area balance (pump)
W.QOHPI = BOHPI * Ayear / 1000000; % Oukhit pump area balance (pump inf)
W.QFJFI = BFJFI * Ayear / 1000000; % Fezna Jorf pump area balance (flood)
W.QFJKI = BFJKI * Ayear / 1000000; % Fezna Jorf Khettara area balance 
                                   % (not the same as pump area)
W.QHKI  = BHKI * Ayear / 1000000;  % Hannabou area balance (khettara)
W.QHFI  = BHFI * Ayear / 1000000;  % Hannabou area balance (flood)
W.QSI   = sum(BSI,2) * Ayear / 1000000;% Sandinfiltration (flood)

W.IN.PInfiltration = mean(W.QFJPI) + mean(W.QOHPI);% [Mm3 / y]
W.IN.FInfiltration = mean(W.QHFI) + mean(W.QFJFI) + mean(W.QSI); % [Mm3 / y]
W.IN.KInfiltration = mean(W.QHKI)+ mean(W.QFJKI) + mean(W.QFKI);  % [Mm3 / y]
W.OUT.QPump        = mean(W.QFJP) + mean(W.QOHP);   % [Mm3 / y]
W.ControlWell = W.TotWell - (W.IN.PInfiltration + W.IN.FInfiltration + ...
    W.IN.KInfiltration + W.OUT.QPump);

%% Evapotransparation
% This part of the model will calculate the ET. It will also look up the ET
% in the areas with infiltration. Sometimes infiltration takes place in
% deeper layers. This is accounted for in the loop, it will look up the
% cells at the top surface, because ET only takes place the top soil (in
% this model at least). ET from the pumping and draining is not taken into
% account. The S.InfFact makes sure 60 - 70 % of the water from pumping and
% draining is not returning to model. This is of course infiltration, but
% this will not be part of the balance for evaporation.

EV = zeros(size(B));    ETFJP = zeros(size(B));
ETFJK = zeros(size(B)); ETHaK = zeros(size(B));
ETOHP = zeros(size(B));

for iE = 1:length(B)
    Eva = strmatchi('ET',B(iE).label);
    EV(iE) = sum(B(iE).term{Eva}(:));
    ET = B(iE).term{Eva};
    
    % FJPI.Idx == FJFI.Idx %
    ZEROSFJPI = zeros(size(gr.AREA3));   % zero matrix all layers
    ZEROSFJPI(pumpAreas(FJPI).Idx) = 1;  % 1 where area is
    ZEROSFJPI = sum(ZEROSFJPI,3);        % sum to 1 layer
    ETFJP(iE) = sum(ET(find(ZEROSFJPI)));% get the ET at top layer
    
    ZEROSOHPI = zeros(size(gr.AREA3));   % zero matrix all layers
    ZEROSOHPI(pumpAreas(OHPI).Idx) = 1;  % 1 where area is
    ZEROSOHPI = sum(ZEROSOHPI,3);        % sum to 1 layer
    ETOHP(iE) = sum(ET(find(ZEROSOHPI)));% get the ET at top layer
    
    ZEROSFJKI = zeros(size(gr.AREA3));   % zero matrix all layers
    ZEROSFJKI(pumpAreas(FJKI).Idx) = 1;  % 1 where area is
    ZEROSFJKI = sum(ZEROSFJKI,3);        % sum to 1 layer
    ETFJK(iE) = sum(ET(find(ZEROSFJKI)));% get the ET at top layer
    
    ZEROSHKI = zeros(size(gr.AREA3));    % zero matrix all layers
    ZEROSHKI(pumpAreas(HKI).Idx) = 1;    % 1 where area is
    ZEROSHKI = sum(ZEROSHKI,3);          % sum to 1 layer
    ETHaK(iE) = sum(ET(find(ZEROSHKI))); % get the ET at top layer
end

W.OUT.ET = mean(EV) * Ayear / 1000000;    % [Mm3 / y]
W.ETFJP  = mean(ETFJP) * Ayear / 1000000; % [Mm3 / y]
W.ETFJK  = mean(ETFJK) * Ayear / 1000000; % [Mm3 / y]
W.ETHaK  = mean(ETHaK) * Ayear / 1000000; % [Mm3 / y]
W.ETOHP  = mean(ETOHP) * Ayear / 1000000; % [Mm3 / y]

%% Results of Wells and ET
% here the results of pumping, infiltration and ET are presented. First of
% the entire model, then of local areas. It will always present the
% results, even if the feature is not selected (it will present 0 values if
% not selected). ET can also be on, but present 0 values. Then there is
% simply no evaporation in this area. There is of course the 60 - 70 %
% evaporation from the pumping and draining, but that is not taken into
% account here.

fprintf('Total pumping, infiltration and ET in the model:\n');
fprintf('Pumping               [Mm3/y] -->'); 
    fprintf('%12.2f',W.OUT.QPump); fprintf('\n');
fprintf('Pumping infiltration  [Mm3/y] -->'); 
    fprintf('%12.2f',W.IN.PInfiltration); fprintf('\n');
fprintf('Khettara infiltration [Mm3/y] -->'); 
    fprintf('%12.2f',W.IN.KInfiltration); fprintf('\n');
fprintf('Flood infiltration    [Mm3/y] -->'); 
    fprintf('%12.2f',W.IN.FInfiltration); fprintf('\n');
fprintf('Evapotransparation    [Mm3/y] -->'); 
    fprintf('%12.2f',W.OUT.ET); fprintf('\n\n');
    
fprintf('Pumping and flood infiltration and ET in Fezna and Jorf (pump area)\n');
fprintf('Pumping extraction   [Mm3/y] -->'); 
    fprintf('%12.2f',mean(W.QFJP)); fprintf('\n');
fprintf('Pumping infiltration [Mm3/y] -->'); 
    fprintf('%12.2f',mean(W.QFJPI)); fprintf('\n');
fprintf('Flood infiltration   [Mm3/y] -->'); 
    fprintf('%12.2f',mean(W.QFJFI)); fprintf('\n');
fprintf('ET                   [Mm3/y] -->'); 
    fprintf('%12.2f',W.ETFJP); fprintf('\n');
fprintf('Khettara Infiltration[Mm3/y] -->'); 
    fprintf('%12.2f',mean(W.QFKI)); fprintf('\n');
fprintf('Balance FJ pump area [Mm3/y] -->');
    fprintf('%12.2f',(mean(W.QFJP) + mean(W.QFJPI)) + mean(W.QFJFI) +...
        W.ETFJP + mean(W.QFKI)); fprintf('\n\n');

fprintf('Khettara infiltration and ET in Fezna and Jorf (khettara area)\n');
fprintf('Khettara infiltration    [Mm3/y] -->'); 
    fprintf('%12.2f',mean(W.QFJKI)); fprintf('\n');
fprintf('ET                       [Mm3/y] -->'); 
    fprintf('%12.2f',W.ETFJK); fprintf('\n');
fprintf('Balance FJ khettara area [Mm3/y] -->');
    fprintf('%12.2f',(mean(W.QFJKI) + W.ETFJK)); fprintf('\n\n');
        
fprintf('Infiltration in Hannabou\n');
fprintf('Khettara infiltration [Mm3/y] -->'); 
    fprintf('%12.2f',mean(W.QHKI)); fprintf('\n');
fprintf('Flood infiltration    [Mm3/y] -->'); 
    fprintf('%12.2f',mean(W.QHFI)); fprintf('\n');
fprintf('ET                    [Mm3/y] -->'); 
    fprintf('%12.2f',W.ETHaK); fprintf('\n');
fprintf('Flood infiltration    [Mm3/y] -->'); 
    fprintf('%12.2f',(mean(W.QHKI) + mean(W.QHFI)) + ...
        W.ETHaK); fprintf('\n\n');
    
fprintf('Pumping, ET and infiltration in Oukhit\n');    
fprintf('Pumping extraction   [Mm3/y] -->'); 
    fprintf('%12.2f',mean(W.QOHP)); fprintf('\n');
fprintf('Pumping infiltration [Mm3/y] -->'); 
    fprintf('%12.2f',mean(W.QOHPI)); fprintf('\n');
fprintf('ET                   [Mm3/y] -->'); 
    fprintf('%12.2f',W.ETOHP); fprintf('\n');
fprintf('Balance Oukhit pump area [Mm3/y] -->');
    fprintf('%12.2f',(mean(W.QOHP) + mean(W.QOHPI)) + W.ETOHP); 
    fprintf('\n\n');   
    
fprintf('Additional river infiltration upstream Fezna\n');  
fprintf('Additional River Infiltration  [Mm3/y] -->'); 
    fprintf('%12.2f',mean(W.QSI)); fprintf('\n');    
    
%% Volume in the aquifer

WaterLvl = H(end).values(:,:,1) - gr.Z(:,:,end);
WaterVol = WaterLvl .* gr.AREA;
WaterVol(isnan(WaterVol)) = 0;
W.WaterVolume = sum(WaterVol(:)) * S.sy / 1000000;
    
%% Waterbalance:

% A table in of the in and out flows (storage included)
IN = struct2table(W.IN);
OUT = struct2table(W.OUT);
display(IN); display(OUT);

% Total in and outflows (without flow out and in of the storage)
TotalIN = W.IN.PInfiltration + W.IN.KInfiltration +...
    W.IN.FInfiltration + W.IN.QNorth + W.IN.QRecharge + W.IN.QRiver;
TotalOUT = W.OUT.QPump + W.OUT.QKhettara + W.OUT.QSouth + W.OUT.ET;

fprintf('Total inflow      [Mm3/y] --> %12.2f\n',TotalIN);
fprintf('Total outflow     [Mm3/y] --> %12.2f\n',TotalOUT);
fprintf('Change in Storage [Mm3/y] --> %12.2f\n',mean(W.StorageSum));
fprintf('                              Positive means water flowing out of the storage\n');
fprintf('unknown           [Mm3/y] --> %12.2f\n',(TotalIN + TotalOUT) + mean(W.StorageSum));
fprintf('                              This value should be zero\n');
fprintf('Total water storage [Mm3] --> %12.2f\n\n',W.WaterVolume);

%% Control of values by checking Budget file
TotS = 0;
TotQ = 0;

for iS = 1:length(B)
    Storage = sum(B(iS).term{1}(:));
    TotS = TotS + Storage;

    Qin = (B(iS).term{2} + B(iS).term{6} + B(iS).term{7} + B(iS).term{8} + B(iS).term{9} + B(iS).term{10});
    Qin = sum(Qin(:));
    TotQ = Qin + TotQ;
end

DeltaStor = TotS / numel(B) * Ayear / 1000000;
DeltaQ    = TotQ / numel(B) * Ayear / 1000000;

if (DeltaStor + DeltaQ) > 0.0001 || (DeltaStor + DeltaQ) < -0.0001
    error('Budget file incorrect. Storage and flows do not match')
elseif (DeltaStor - mean(W.StorageSum)) > 0.0001 || ...
        (DeltaStor - mean(W.StorageSum)) < -0.0001
    error('Storage from the budget file and the waterbalance do not match')
elseif (TotalIN + TotalOUT) + mean(W.StorageSum) > 0.01 || ...
        (TotalIN + TotalOUT) + mean(W.StorageSum) < -0.01
    error('Unknown balances in the waterbalance.')
elseif W.ControlWell > 0.0001 || W.ControlWell < -0.0001
    error('Balance in budget file well is incorrect.')
elseif ~S.Floods
    if abs(W.IN.FInfiltration) > 0.00001 ||  abs(W.IN.FInfiltration) < -0.00001
        error('There is still flood inf even though S.Floods is not true');
    end
elseif ~S.Pumping
    if abs(W.OUT.QPump) > 0.00001 ||  abs(W.OUT.QPump) < -0.00001
        error('There is still pumping even though S.pumping is not true');
    elseif abs(W.IN.PInfiltration) > 0.00001 ||  ...
            abs(W.IN.PInfiltration) < -0.00001
        error('There is still pumping inf even though S.pumping is not true');
    end
elseif ~S.Khettaras
    if abs(W.OUT.QKhettara) > 0.00001 ||  abs(W.OUT.QKhettara) < -0.00001
        error('There is still Khet discharge even though S.Khettaras is not true');
    elseif abs(W.IN.KInfiltration) > 0.001 ||  ...
            abs(W.IN.KInfiltration) < -0.01
        error('There is still khet inf even though S.Khettaras is not true');
    end
end

%% Save the waterbalance, so it can be used in a next run.
save Waterbalance W

%% Plot the khettaras length at specified time
if S.Khettaras && S.KhetPlot
    % plot the khettara length of the time step and the part of khettara
    % which is draining water.
    
    axProps = {'nextplot','add','xGrid','on','yGrid','on'}; % default figure porporties
    formationColors = 'ymcg';
    XLIM = [3.52*10^5 gr.xGr(end)];
    YLIM = [gr.yGr(end) 3.492*10^6];
    
    figure('name','head relative ground surface','pos',screenPos(0.75),'renderer','zbuffer');
    axes(axProps{:},'xlim',XLIM,'ylim',YLIM);
    xlabel('x UTM [m]'); ylabel('y UTM [m]');
    
    try
        title(sprintf('Head relative to ground surface, with khettaras L = %d m in the year %d ',S.Radius,S.Years));
    catch
        title('head relative to ground surface, with khettaras length');
    end

    [c,h] = gr.contourf(H(end).values(:,:,1),780:1:850);  clabel(c,h,'color','y'); %-gr.Z(:,:,1)

    % Add white contour showing the zero line where head is at ground surface
    gr.contour (H(end).values(:,:,1)-gr.Z(:,:,1) ,[0 0],'w','lineWidth',1);

    hb = colorbar; set(get(hb,'title'),'string','\Phi - DEM [m]');
    set(gcf, 'renderer', 'zbuffer');
    
    [KH] = plot(gr.XM([khettaras.Idx]),gr.YM([khettaras.Idx]),'ko','markersize',3,'linewidth',3);
    [KHW]= plot(gr.XM(IdxKhet),gr.YM(IdxKhet),'mo','markersize',3,'linewidth',3);
    legend([KH KHW],'Khettara','Khettara water extraction');
end

%% Used to check the results
if S.Control
    % Here you can plot different areas, to see where they are located in
    % the model. this will provide usefull feedback to see whether the
    % areas are correctly located.
    
    figure; hold on;
    [ibound] = plot(gr.XM(IBOUND==0),gr.YM(IBOUND==0),'yo');
    idx = pumpAreas(FJP).Idx(std(WFJP,1,2) > 0.0001);
    intersect(find(IBOUND~=0),pumpAreas(FJP).Idx);
    [Contourareas] = plot(pumpAreas(2).vertex(:,1),pumpAreas(2).vertex(:,2)); 
    [PumpA] = plot(gr.XM(intersect(find(IBOUND~=0),pumpAreas(FJP).Idx)),gr.YM(intersect(find(IBOUND~=0),pumpAreas(FJP).Idx)),'ro');
    [PumpEmpty] = plot(gr.XM(idx),gr.YM(idx),'mo');
    [OutCrop] = plot(gr.XM(vertcat(outcrops.Idx)),gr.YM(vertcat(outcrops.Idx)),'ko');
    % [cc,hh] = contour(gr.Xm,gr.Ym,H(1).values(:,:,4) - gr.Z(:,:,end),0:20);
    % clabel(cc,hh);
    
    for i=1:length(khettaras)
        if any(intersect(pumpAreas(FJKI).Idx,khettaras(i).Idx))
            [KhetOutIbound] = plot(khettaras(i).vertex(:,3),khettaras(i).vertex(:,4),'ro');
            khettaras(i).name
        else
            [KhetInIbound] = plot(khettaras(i).vertex(:,3),khettaras(i).vertex(:,4),'go');
        end
    end
    
    plot(gr.XM(pumpAreas(HKI).Idx),gr.YM(pumpAreas(HKI).Idx),'ko')
    legend([ibound Contourareas PumpA PumpEmpty OutCrop KhetOutIbound...
        KhetInIbound],'IBOUND = 0','Contour pump area','Cells in pump area',...
        'Pump area with dry aquifer','Outcrops','khettaras draining out of model',...
        'Khettaras in Model');
    
    
    % Here you can check whether the total storage is equal to the total
    % outflows (-inflows) in the model.
    TotS = 0;
    TotQ = 0;

    for iS = 1:length(B)
        Storage = sum(B(iS).term{1}(:));
        TotS = TotS + Storage;

        Qin = (B(iS).term{2} + B(iS).term{6} + B(iS).term{7} + B(iS).term{8} + B(iS).term{9} + B(iS).term{10});
        Qin = sum(Qin(:));
        TotQ = Qin + TotQ;
    end
    
    %% control of the khettaras
    
    for i=1:length(khettaras)
        Startx(i) = khettaras(i).P(end).x(2); 
        Starty(i) = khettaras(i).P(end).y(2);
        Eindx(i) = khettaras(i).P(1).x(2);
        Eindy(i) = khettaras(i).P(1).y(2);
        drainH(i)= khettaras(i).P(1).zm;

        EindIdx  = khettaras(i).P(1).idx;
        StartIdx = khettaras(i).P(end).idx;

        IX = khettaras(i).P(1).ix;
        IY = khettaras(i).P(1).iy;

        LAY1(i) = gr.Z(IY,IX,1);
        LAY2(i) = gr.Z(IY,IX,2);
        LAY3(i) = gr.Z(IY,IX,3);
        LAY4(i) = gr.Z(IY,IX,4);
        LAY5(i) = gr.Z(IY,IX,5);
    end

    khet = [Eindy' Eindx' LAY1' LAY2' LAY3' LAY4' LAY5' drainH'];
    khet = array2table(khet);
    khet = sortrows(khet,1,'descend');
    % khet.Properties.VariableNames ={'EindY','EindX','TopElevation','DownElevation','DrainHoogtebegininvlakte'};
    display(khet);

    figure; hold on;
    plot(table2array(khet(:,1)),table2array(khet(:,8)),'bo');
    plot(table2array(khet(:,1)),table2array(khet(:,3)),'r');
    plot(table2array(khet(:,1)),table2array(khet(:,3)),'m');
    plot(table2array(khet(:,1)),table2array(khet(:,4)),'c');
    plot(table2array(khet(:,1)),table2array(khet(:,5)),'y');
    plot(table2array(khet(:,1)),table2array(khet(:,6)),'k');
    plot(table2array(khet(:,1)),table2array(khet(:,7)),'r');
end
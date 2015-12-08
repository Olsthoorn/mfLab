function Run_Model(S)
% This function will run the model assosiated with Koen_scenarios.m.

if S.KhetLength
    % If khetlength = true, this part will model the khettara length over
    % time
    S.KhetInfo = 0;
    S.Radius = 0;
    S.Years = S.StartYKhet; 
    save Scenarios S
    
    for iK = 1:numel(S.StartYKhet:11:S.EindYKhet)
        if S.Random
            S.PercentFlooding = mean(S.y(:,iK));
            S.PercentRecharge = S.y(:,iK);
            S.PercentRiver    = mean(S.y(:,iK));
        end
        
        save Scenarios S
        save NumberKhet iK
        %load NumberOfRuns
        
        for iWa = 1:S.NumberOfKhetRuns
            save Numbercount iWa % this is neccesary because mf_setup clears variables
            mf_setup;
            load Numbercount
            Koen_waterbalance;
            load Waterbalance

            if abs(mean(W.StorageP + W.StorageN)) < S.MaxStorCh && ~S.Random
                break;
            end
        end
        
        load NumberKhet
        load Scenarios

        if iK==1, S.KhetInfo = NaN(numel(S.StartYKhet:11:S.EindYKhet), 15); end
        
        S.KhetInfo(iK,:) = [S.Years -mean(W.QFJP) -mean(W.QOHP) -W.OUT.QKhettara ... 
            -W.OUT.ET W.IN.QRecharge W.IN.QRiver W.IN.FInfiltration ...
            mean(W.StorageSum) W.WaterVolume W.IN.PInfiltration ...
            W.IN.KInfiltration W.IN.QNorth -W.OUT.QSouth S.Radius];
        S.Years = S.Years + S.YearKhet;
        S.Radius = S.Radius + S.Length;
        save Scenarios S
    end
    mf_analyze;
    
    % Process information of the khettara length
    load Scenarios
    Plottable(S.KhetInfo,1,1);
    
elseif S.TimePump
    % This part will simulate pumping rate over time
    if S.Oukhit
        S.PumpVolOuk = 0;
        S.PumpInfo = 0;
        S.Years = S.StartYPump;
        S.PumpVolMax = S.PumpOukhitMax;
        S.PumpVolCh = S.PumpVolChOuk;
        save Scenarios S
    else
        S.PumpVol = 0;
        S.PumpInfo = 0;
        S.Years = S.StartYPump;
        save Scenarios S
    end
    
    for iP = 1:numel(1:S.PumpVolCh:S.PumpVolMax)
        save NumberPump iP
        
        for iWa = 1:S.NumberOfRuns
            save Numbercount iWa % this is neccesary because mf_setup clears variables
            mf_setup;
            load Numbercount
            Koen_waterbalance;
            load Waterbalance

            if abs(mean(W.StorageP + W.StorageN)) < S.MaxStorCh && ~S.Random
                break;
            end
        end
        
        load Scenarios
        load NumberPump

        S.PumpInfo(iP,:) = [S.Years -mean(W.QFJP) -mean(W.QOHP) -W.OUT.QKhettara ... 
            -W.OUT.ET W.IN.QRecharge W.IN.QRiver W.IN.FInfiltration ...
            mean(W.StorageSum) W.WaterVolume W.IN.PInfiltration ...
            W.IN.KInfiltration W.IN.QNorth -W.OUT.QSouth];
        S.Years = S.Years + S.YearPump;
        
        if S.Oukhit
            S.PumpVolOuk = S.PumpVolOuk + S.PumpVolCh;
        else
            S.PumpVol = S.PumpVol + S.PumpVolCh;
        end           
        save Scenarios S
    end
    mf_analyze;    
    % Process information of Pumping
    load Scenarios
    Plottable(S.PumpInfo,2,1);
    
else
    for iWa = 1:S.NumberOfRuns
        if S.Random
            S.PercentFlooding = S.y(iWa);
            S.PercentRecharge = S.y(iWa);
            S.PercentRiver    = S.y(iWa);
            save Scenarios S
        end
        
        fprintf('.. Model run number %d\n',iWa);

        save Numbercount iWa % this is neccesary because mf_setup clears variables
        mf_setup;
        load Numbercount
        Koen_waterbalance;

        load Waterbalance
        load Scenarios
        
        S.ModelInfo(iWa,:) = [iWa -mean(W.QFJP) -mean(W.QOHP) -W.OUT.QKhettara ... 
            -W.OUT.ET W.IN.QRecharge W.IN.QRiver W.IN.FInfiltration ...
            mean(W.StorageSum) W.WaterVolume W.IN.PInfiltration ...
            W.IN.KInfiltration W.IN.QNorth -W.OUT.QSouth];

        if abs(mean(W.StorageP + W.StorageN)) < S.MaxStorCh && ~S.Random
            fprintf('\n\n\n.. Ended at model run %d\n\n\n\n',iWa);
            break;
        end
        save Scenarios S
    end
    mf_analyze;
    Plottable(S.ModelInfo,3,1);
end
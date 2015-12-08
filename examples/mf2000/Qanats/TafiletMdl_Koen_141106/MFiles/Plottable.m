function Plottable(o,a,b)
% will create tables from info of the model and make plots of the results.
% 
% KG 141028

if b
    o = array2table(o);
end

switch a
    case 1
        o.Properties.VariableNames = {'Run','PumpingFeznaJorf',...
            'PumpingOukhit','KhettaraDischarge','ET','Recharge',...
            'RiverInfiltration','FloodInfiltration','ChangeStorage',...
            'StorageVolume','PInfiltration','KInfiltration',...
            'BoundaryNorth','BoundarySouth','Length'};
        o.Properties.VariableUnits = {'[T]','[Mm3/y]',...
            '[Mm3/y]', '[Mm3/y]','[Mm3/y]','[Mm3/y]', ...
            '[Mm3/y]','[Mm3/y]','[Mm3/y]',...
            '[Mm3]','[Mm3/y]','[Mm3/y]',...
            '[Mm3/y]','[Mm3/y]','[m]'};
        display(o);
        if b
            save KhettaraInformation o
        end
        Length = table2array(o(:,15));
    case 2
        o.Properties.VariableNames = {'Run','PumpingFeznaJorf',...
            'PumpingOukhit','KhettaraDischarge','ET','Recharge',...
            'RiverInfiltration','FloodInfiltration','ChangeStorage',...
            'StorageVolume','PInfiltration','KInfiltration',...
            'BoundaryNorth','BoundarySouth'};
        o.Properties.VariableUnits = {'[T]','[Mm3/y]','[Mm3/y]','[Mm3/y]', ...
            '[Mm3/y]','[Mm3/y]','[Mm3/y]','[Mm3/y]','[Mm3/y]','[Mm3]',...
            '[Mm3/y]','[Mm3/y]','[Mm3/y]','[Mm3/y]'};
        display(o);
        if b
            save PumpingInformation o
        end
    case 3
        o.Properties.VariableNames = {'Run','PumpingFeznaJorf',...
            'PumpingOukhit','KhettaraDischarge','ET','Recharge',...
            'RiverInfiltration','FloodInfiltration','ChangeStorage',...
            'StorageVolume','PInfiltration','KInfiltration',...
            'BoundaryNorth','BoundarySouth'};
        o.Properties.VariableUnits = {'[T]','[Mm3/y]','[Mm3/y]','[Mm3/y]', ...
            '[Mm3/y]','[Mm3/y]','[Mm3/y]','[Mm3/y]','[Mm3/y]','[Mm3]',...
            '[Mm3/y]','[Mm3/y]','[Mm3/y]','[Mm3/y]'};
        display(o);
        if b
            save ModelInformation o
        end
end

Year        = table2array(o(:,1));
PumpFJ      = table2array(o(:,2));
PumpOH      = table2array(o(:,3));
DischargeK  = table2array(o(:,4));
ET          = table2array(o(:,5));
Recharge    = table2array(o(:,6));
RiverInf    = table2array(o(:,7));
FloodInf    = table2array(o(:,8));
ChangeStor  = table2array(o(:,9));
StorVol     = table2array(o(:,10));
PInf        = table2array(o(:,11));
KInf        = table2array(o(:,12));
NorthB      = table2array(o(:,13));
SouthB      = table2array(o(:,14));

if a == 1
    figure; hold on; plot(Year,DischargeK); title('Khettara Discharge');
    xlabel('Year'); ylabel('Discharge [Mm3/y]');
    figure; hold on; plot(Length,DischargeK); title('Khettara Discharge');
    plot(Length,KInf,'r'); legend('Discharge','Infiltration');
    xlabel('Khettara length'); ylabel('Discharge [Mm3/y]');
else
    figure; hold on; plot(Year,DischargeK,'r'); plot(Year,PumpFJ,'b');
    plot(Year,PumpOH,'m'); xlabel('Year'); ylabel('Discharge [Mm3/y]');
    plot(Year,PInf,'g'); plot(Year,KInf,'k'); 
    title('Pumping and Khettara Discharge');
    legend('Khettara Discharge','Pumping Discharge Fezna Jorf',...
        'Pumping Discharge Oukhit','Infiltration');
end

figure; hold on;
title('Discharges in the area');
    plot(Year,ET,'c');
    plot(Year,Recharge,'c--');
    plot(Year,RiverInf,'m');
    plot(Year,FloodInf,'m--');
    plot(Year,DischargeK,'r');
    plot(Year,PumpFJ,'b');
    plot(Year,PumpOH,'g');
    plot(Year,PInf,'b--');
    plot(Year,KInf,'r--');
    plot(Year,NorthB,'k');
    plot(Year,SouthB,'k--');
    legend('ET','Recharge','River Infiltration','Flood Infiltration',...
        'Discharge Khettara','Pumping Discharge Fezna Jorf',...
        'Pumping Discharge Oukhit','Pump infiltration',...
        'Khettara infiltration','North Boundary','South Boundary');
    xlabel('Year'); ylabel('Discharge [Mm3/y]');

figure;
    subplot(2,2,1);
    title('Pumping + khettara discharge');
        if a == 1
            plot(Year,DischargeK);
        else
            hold on; plot(Year,DischargeK,'r');  plot(Year,PumpFJ,'b');
            plot(Year,PumpOH,'m'); xlabel('Year'); 
            ylabel('Discharge [Mm3/y]');
            legend('Khettara Discharge','Pumping Discharge Fezna Jorf',...
                'Pumping Discharge Oukhit');
        end
        xlabel('Year'); ylabel('Discharge [Mm3/y]');
    subplot(2,2,2); plot(Year,StorVol);
        title('Volume in storage');
        xlabel('Year'); ylabel('Storage Volume [Mm3]');
    subplot(2,2,3); plot(Year,ChangeStor);
        title('Change in storage');
        xlabel('Year'); ylabel('Change in Storage [Mm3/y]');
    subplot(2,2,4); plot(Year,Recharge); 
        title('Recharge in the area');
        xlabel('Year'); ylabel('Recharge [Mm3/y]');
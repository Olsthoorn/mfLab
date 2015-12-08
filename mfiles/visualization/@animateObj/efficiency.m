function well = efficiency(o,well,iCompTEmp)
    %ANIMATEOBJ/EFFICIENCY -- compute and show the efficiency of all wells
    %versus time
    %
    % USAGE  well = animateObj.efficiency(well,iCompTemp)
    %  iCompTemp is the compoment number associated with temperature
    %  well i of class wellObj
    %  The computed efficiency is added to the well in the output.
    
    % Get output temperatures for the wells
    well = well.setCout(o.UCN{iCompTEmp},iCompTEmp);

    % Then compute the efficiency of the wells and add this info to their
    % UserDta.
    well = well.efficiency(iCompTEmp);

    % Then show this efficiency graphically
    figure('Name','Temperature efficiency of the wells');
    axes('nextplot','add','xgrid','on','ygrid','on');
    if well(end).t(end)>3650,
        xlabel('time [y]');
        tfac = 365.24;
    else
        xlabel('time [d]');
        tfac = 1;
    end
    ylabel('temp. efficiency');
    title('Output temperature efficiency of the wells');
    
    % For each well draw the line
    leg{numel(well),1}='';
    type ='WC';  % warm or cold
    % UserData.type == 'D' for doublet and 'M' for monowell.
    for iw=1:numel(well)
        if ~isfield(well(iw).UserData,'type'), well(iw).UserData.type = type; end
        if well(iw).UserData.type   == 'D',   lineWidth = 2;    else lineWidth =   1; end
        if well(iw).UserData.isCold == false, lineStyle = '--'; else lineStyle = '-'; end
            leg{iw}=sprintf('well%02d %s %s scr %.0f %.0f',...
                well(iw).nr, ...
                well(iw).UserData.type, ...
                type(well(iw).UserData.isCold+1), ...
                well(iw).z([1 end]));
            plot(well(iw).t(2:end)/tfac,well(iw).UserData.efficiency,well(iw).UserData.color,...
                'lineStyle',lineStyle,'lineWidth',lineWidth)
    end
    legend(leg,4);
end
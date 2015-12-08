function o = efficiency(o,iComp)
    %WELLOBJ/EFFICIENCY -- computes and plots seasonal thermal efficiency
    % and ads the result to the  UserData field of the wellObj.
    %
    % well = well.efficiency(iComp)
    % iComp must be the species' index of the temperature in the simulation
    % and well must have the temperature on board in well(iw).Cout(iComp,:)
    % well = well.setCout(T,iComp) must, therefore, have already been
    % excecuted
    %
    % TO 120913
    
    aYear = 365.24;
    aWeek = 7;       % space around a year
    
    if isempty(o(1).Cout) || all(isnan(o(1).Cout(iComp,:)))
        error('%s: Cout not set, apply well=well.setCout(C,iComp) first',mfilename);
    end
    
    % Plot output concentration of wells versus time  
    % Because we compute the thermal efficiency on a well by well basis, we
    % cannot use the information of other wells. Hence, we compute the
    % efficiency by comparing the output of a well with the injecion one
    % season ago. But to make this computation continuous and its
    % implementation and interpretation more straightforward, we compare
    % for any whole year, irresepctive of its start, the total extracted
    % energy with the total injected energy in the same well. That is we
    % compare 6 month injection during a year with 6 month extraction
    % during the same year, by intgration over that year.
    for iw=1:numel(o)
        extraction = o(iw).Q<0;
        injection  = ~extraction;
        Eout =    -extraction .* o(iw).Q.* o(iw).Dt .* (o(iw).Cout(iComp,:)-o(iw).UserData.taverage);
        Ein  =     injection  .* o(iw).Q.* o(iw).Dt .* (o(iw).C(   iComp,:)-o(iw).UserData.taverage);
        
        % Get distance between times that are a seaon apart, i.e. 6 months
        % apart
        % notice that we have stress periods of one month.
        iEnd   = 1:numel(o(iw).Dt);
        iStart = NaN(size(iEnd));
        for it=iEnd
            % Istrt is the start month, which must be 1 year before the end
            % month. 
            istrt = find(o(iw).t>=o(iw).t(it)-aYear-aWeek & o(iw).t<o(iw).t(it)-aYear+aWeek);
            if ~isempty(istrt),
                % The start month number is 1 more than the
                iStart(it)=istrt+1;
            else
                iStart(it)=NaN;
            end
        end        
        if all(isnan(iStart)) % Then we choose an iEnd shorter than a year
            msgId = 'mfLab:seasonalATESwellEfficiency:simPeriodSmallerThanAYear';
            warning('on',msgId);
            beep;
            warning(msgId,...
                ['Time series for wells must be > 1 year to be able to compute its seasonal efficiency!\n',...
                'Remedy: make sure the total simulation time in the PER sheet is > 365 days, rather > 730 d']);
            warning('off',msgId);
        end        
        o(iw).UserData.efficiency = NaN(size(o(iw).Dt));
        for it=1:numel(iStart)
            if ~isnan(iStart(it))
                % integrate Eout/ ingegrate Ein over the year startin at
                % iStart(it)
                o(iw).UserData.efficiency(it)=sum(Eout(iStart(it):iEnd(it)))./sum(Ein(iStart(it):iEnd(it)));
                if o(iw).UserData.efficiency(end)<0
%                    error('well(%d) efficiency negative !',iw);
                end
            end
        end
    end
end

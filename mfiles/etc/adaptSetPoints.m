function well=adaptSetPoints(well,tTemp,nYr,~)
%ADAPTSETPOINTS  adapt set points for ATES systems
%    ATES stands for "Aquifer thermal energy storate"
%
%USAGE:
%    well=adaptSetPoints(well,tTemp,nYr,~)
%
% See also: mflab/examples/swt_v4/ATES-WKO/Utrecht-NL/Utrecht5YrSetpUpdate/mf_adapt
%
% Adapt setting points of wells for cooling and heating such that the
% storage and use of the subsurface is in balance over a 5-year period.
%
% Function is specific to the ATES problem of Utrecht as worked out
% under examples mflab/examples/swt_v4/ATES-WKO/Utrecht-NL/
%
% Temp is the [time temp] vector of the air.
%
%PROCEDURE:
%   Each well has 4 temperature setpoints, two of which are used if it is a warm
%   well and two in case it is a cold well.
%   The temperature setpoints are as follows
%       Tqmax_c (tqMc) = temperature at which Qcooling is at its maximum
%       Tqmin_c (tqmc) = temperature at which Qcooling is miniumum (switches on)
%       Tqmax_h (tqMh) = temperature at which Qheating is at its maximum
%       Tqmin_h (tqmh) = temperature at which Qheating is at its minimum
%
%   tqmc and tqmh are optimized such that total heating demand and cooling
%   demand match over the considered number of years.
%
%   Tqmax_c and Tqmax_h are given and remain fixed
%   Tqmin_c and Tqmin_h are reset.
%
%   nYr is the number of years over which the setpoints are computed.

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

verbose = nargin>3;

if nargin<4 || isempty(nYr), nYr=5; end % default

f=well.UserData.QfracMin;   % minimum well flow as fraction of its capacity

%% Get years in time series
DV=datevec(tTemp(:,1)); Y=unique(DV(:,1));

%% Allocate vector to store the swicht temperatures
Tqmin_c = NaN(length(Y)-nYr,1);  % setpoints for Tqmin_c
Tqmin_h = NaN(length(Y)-nYr,1);  % setpoints for Tqmin_h

%% Compute length(Y)-nYr yearly temperature setpoints.
ii=0;
for iy=(nYr+1):length(Y)
    
    %% get temperature for the current nYr period
    temp=tTemp(DV(:,1)>=Y(iy-nYr) & DV(:,1)<Y(iy),2);
    
    %% COOLING switching on setpoint

    T0=0;  % initial value of tqmc2. Must be lower than anypossible temperature.
    
    tqMc=well.UserData.Tqmax_c; % temp at which cooling from ATES is at its maximum (fixed)
    
    %% The optimal temperature will be somewhere between tqmc1 and tqmc2
    tqmc1=tqMc;  % top   trial temperature at which cooling is switched on
    tqmc2=T0;    % lower trial temperature at which cooling is swithced on
    tqmc =0.5*(tqmc1+tqmc2); % initial trial temp, midway between tqmc1 and qqmc2
    
    %% Use the regula falsi method to compute tqmc, a scalar which we put in a
    %  vector with year values upon breaking out of the while loop
    for k=1:100
        
        % Average daily and yearly extraction for cooling over the considered period 
        QdCooling=well.UserData.QdMax*(temp>tqmc).*min(1,f+(1-f)*(temp-tqmc)/(tqMc-tqmc));
        
        QyCooling=sum(QdCooling)/nYr;  % average QyCooling depends on setpoint
                                       % and must equal predescribed
                                       % well.Qy/2
        
        if QyCooling<well.UserData.Qy % If total QyCooling < desired total extraction for cooling
            tqmc1=tqmc;               % reduce temp at which cooling starts (more cooling)
            tqmc=0.5*(tqmc1+tqmc2);
        else                          % if total QyCooling > desired
            tqmc2=tqmc;               % raise temp at which cooling starts (less cooling days)
            tqmc =0.5*(tqmc1+tqmc2);
        end
        %% The setpoint is such that the QyCooling in the nYr years previous
        %  to the current year equals Qy/2.

        if tqmc1-tqmc2<0.000001;  % Break off criterion
            % this checks if the flow is correct.
            % This must be done every nYr period
            if verbose
                fprintf('well %d %10s year = %4d, QCyr=%12.0f QyCooling=%12.0f tqmc=%.2f\n',...
                    well.nr,well.name,iy,well.UserData.Qy,QyCooling,tqmc);
            end
            Tqmin_c(iy-nYr)=tqmc; % one value per year
            break;
        end
    end
    
    if verbose
        if k==100, s='not'; else s=''; end
        fprintf('%s: Loop setting tqmc %d finished, k=%d\n',mfilename,s,k)
    else
        ii=ii+1;
        if rem(ii,50)  ==0, fprintf('.'); end
    end
    
    %% HEATING switch-on setpoint
    
    T0=20;             % Initial temp at which heating is switched on (heating is minimum)
    tqMh = well.UserData.Tqmax_h; % Temp at which heating is maximum (fixed)
    
    %% The optimal temperature is somewhere between tqmh1 and tqmh2
    tqmh1=tqMh;              % upper trial temperature (set equal to max heating temp)
    tqmh2=T0;                % lower trial temperature (set equal to min heating temp initial)
    tqmh =0.5*(tqmh1+tqmh2); % start in the middle and adapt iteratively
    
    %% Compute setpoint by method of regula falsi
    for k=1:100
        QdHeating=well.UserData.QdMax*(temp<tqmh).*min(1,f+(1-f)*(temp-tqmh)/(tqMh-tqmh));
        QyHeating=sum(QdHeating)/nYr;
        
        if QyHeating<well.UserData.Qy  % too little flow compared to design value
            tqmh1=tqmh;                  % raise temp at which heating starts to get more heating days in a year 
            tqmh=0.5*(tqmh1+tqmh2);
        else % Too much heating flow compared to desired design value
            tqmh2=tqmh;    % lower temp at which heating starts to get less heating days in a year
            tqmh =0.5*(tqmh1+tqmh2);
        end

        if abs(tqmh1-tqmh2)<0.000001; % break off criterion
            % this checks if flows are correct and must be done per period
            if verbose
                fprintf('well %d %10s year = %4d, QHyr=%12.0f QyHeating=%12.0f tqmh=%.2f\n',...
                    well.nr,well.name,iy,well.UserData.Qy,QyHeating,tqmh);
            end
            Tqmin_h(iy-nYr)=tqmh; % one value per year
            break;
        end
    end
    if verbose
        if k==100, s='not'; else s=''; end
        fprintf('%s: Loop setting tqmh %d finished, k=%d\n',mfilename,s,k)
    else
        ii=ii+1;
        if rem(ii,  50)==0, fprintf('.');       end
    end
end

fprintf('%d\n',ii);
%% Finish off
well.UserData.Tqmin_c=[Y((nYr+1):end) Tqmin_c];      % join year with setpoint value
well.UserData.Tqmin_h=[Y((nYr+1):end) Tqmin_h];      % join year with setpoint value



 
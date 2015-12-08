function Par = parChange(Par,changeBy)
%%PARCHANGE -- multiply the paramters with the values in changeBy for as far
%              as Par.use is true or nonzero
%
% This function changes the parameter values in Par by multiplying
% the values by those in changeBy.
% prepares for calibration a model by non-linear optimization
% routines of Matlalb, such as p = lsqnonlin(@fun,p0);
%
% USAGE:
%   Par = parChange(Par,changeBy)
%    Par is set of all parameters, cell array with values
%       {NAME   VALUE LB UB LN USE; ...}
%        NAME is name of parameter
%        VALUE is its value
%        LN    flag indicating to use ln-transformation
%        USE = flag indicating to calibrate this parameter
%        for flags use eithe zero or nonzer or true or false, don'tmix
%   changeBy is the same, with or without LB and UB or LN
%        only its name will be checked against that of par
%        and its value will be used to change Par if use is on
%        if ln is on  Par is multiplied by changeBy
%        if ln is off changeBy is added to Par
%
%    Par should be pase
% TO 130620

    parCheck(Par,'Par')
    parCheck(changeBy,'changeBy');

    for ip=1:size(Par,1)
        if (Par{ip,end})  % used
            if ~strcmp(Par{ip,1},changeBy{ip,1})
                error('Par{%d,1}~=changeBy{%d,2} <<%s~=%s>>',...
                    ip,ip,Par{ip,1},changeBy{ip,2});
            end
            if Par{ip,end-1}
                Par{ip,2} = Par{ip,2}*changeBy{ip,2};
            else
                Par{ip,2} = Par{ip,2}+changeBy{ip,2};
            end        
        else
            % par is lef as it was
        end
    end
end

function parCheck(Par,parname)
    if ~iscell(Par) ...
            || ~all( cellfun(@ischar,Par(:,1))) ...
            || ~all( cellfun(@isscalar,Par(:,2))) ...
            || ~ (all( cellfun(@isnumeric,Par(:,end-1))) || ...
                  all( cellfun(@islogical,Par(:,end-1))))   ...
            || ~ (all( cellfun(@isnumeric,Par(:,end  ))) || ...
                  all( cellfun(@islogical,Par(:,end  ))))
        error(['%s must be cell array with 3 collumns {VARNAME VALUE LB UB LN USE; ...}\n',...
               'with VARNAME of type  char, VALUE a scalar and LN/USE of type numeric or logical'],...
               parname);
    end
end
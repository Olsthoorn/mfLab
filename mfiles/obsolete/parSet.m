function [P,lb,ub,pars] = parSet(Par)
%%PARSET -- get P to be calibrated and prepare Par for calibration
%
% This function prepares for calibration a model by non-linear optimization
% routines of Matlalb, such as p = lsqnonlin(@fun,p0);
% Make sure that Par is global in fun and the workspace, so that all
% parameters are known to fun. In fun reassemble the set using also the
% parameters passed by lsqnonlin.

% USAGE:
%   [P,lb,ub,par] = parSet(Par)
%    P = list of ln transformed parameters to be calibrated
%    lb= lower bound of transformed parameter
%    ub= upper bound of transformed parameter
%    Par is set of all parameters, cell array with values
%       {name   value use; ...}
%        name is name of parameter
%        value is its value
%        use = true of used in calibation
%        use = false if used by model but not calibrated
%
%    Par should be pase
% TO 130620

if ~iscell(Par) ...
        || ~all( cellfun(@ischar,Par(:,1))) ...
        || ~all( cellfun(@isscalar,Par(:,2))) ...
        || ~ (all( cellfun(@isnumeric,Par(:,end))) || ...
              all( cellfun(@islogical,Par(:,end))) ...
             )
    error(['Par must be cell array with 3 collumns {varName, value, use; ...}\n',...
           'with varName of type char, value a scalar and use of type numeric or logical']);
end
    
% get parameters
P = NaN(size(Par,1),1);
lb= NaN(size(Par,1),1);
ub= NaN(size(Par,1),1);

pars{1,2*size(Par,1)} = '';

for ip=1:size(Par,1)
    pars{2*ip-1} = Par{ip,1};
    pars{2*ip  } = Par{ip,2};
    if Par{ip,end} % if parameter needs to be calibrated
        if Par{ip,end-1}
            % multiplyer, log par, initially 1
            P(ip)  = 0;     % i.e. log(Par{ip,2}/Par{ip,2}); 
            lb(ip) = log(Par{ip,3}/Par{ip,2});
            ub(ip) = log(Par{ip,4}/Par{ip,2});
        else % linear parameter
            P(ip)  = Par{ip,2};   % multiplyer, log par, initially 1
            lb(ip) = Par{ip,3};
            ub(ip) = Par{ip,4};
        end
    end
end
P( isnan(P))=[];
lb(isnan(P))=[];
ub(isnan(P))=[];


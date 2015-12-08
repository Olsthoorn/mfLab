function [pars,Par] = parNew(Par,P)
%%PARNEW -- get replacing P in Par by calibrated values P
% works in conjunction wit parSet
%
% This uses ParOld as prepared with parSet and replaces the values in ParOld
% with the calibrated values P. The calibrated values are the ln transformed
% calibrated values, so the replacement is exp(P);
% Works in conjunction with opmization routines of Matlalb,
% such as p = lsqnonlin(@fun,p0);

% USAGE:
%   ParNew = parSet(ParOld,P)
%    P = list of calibrated ln transformed parameters
%        to be replaced in ParOld to form ParNew.
%    Par is the set of all parameters, as cell array
%       with values
%       {name   value use; ...}
%        name is name of parameter
%        value is its value
%        use = true of used in calibation
%        use = false if used by model but not calibrated
%
% TO 130620

if ~iscell(Par) || ~all(isscalar(P))
    error(['Par must be cell array with 3 collumns {varName, value, use; ...}\n',...
           'with varName of type char, value a scalar and use of type numeric or logical.\n',...
           'P must be vector of calibrated values, all scalars']);
end
    
nVar = size(Par,1);
pars = cell(1,size(Par,1));

j=0;
for iv=1:nVar
    if Par{iv,end}
        j=j+1;
        Par{iv,2}=Par{iv,2}*exp(P(j)); % log par
    end
    pars{2*iv-1} = Par{iv,1};
    pars{2*iv  } = Par{iv,2};
end

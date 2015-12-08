function Pars = parShow(cases,varargin)
%%PARSHOW -- shows parameters
% works in conjunction wit parSet, parNew
%
% parameters as cell arrays as used in parSet will be printed
% in the form of {varName use par par par; ...}
% The varnames and use are obtained from the first par given in vararin.
% The other parameters are obtained form subsequent varargin entries.
% Varnames must be the same in the varargins.
%
%
% USAGE:
%   Pars = parShowet(cases Par Par Par Par ...)
%    cases is a cell array of strings denoting the bases pertaning
%    to the Par(1), Par(2) etc.
%    if omitted kind is treated as Par and 'case n' will be used instead
%    as headers of the output columns.
%    Par is the set of all parameters, as cell array
%       with values
%       {name   value use; ...}
%        name is name of parameter
%        value is its value
%        use = true of used in calibation
%        use = false if used by model but not calibrated
%
% TO 130620

if ~iscell(cases)
    error('kind must be a set of n strings where n is length(varargin)');
end

if ~isvector(cases) || ~all(cellfun(@ischar,cases))
    cases=[];
    varargin = [kind varargin];
end

stdP    = varargin{end-1};
uncertP = varargin{end};
varargin(end-1:end) = [];  % varargin only contains the paramter cell arrays

nVar = size(varargin{1},1);
nPar = length(varargin);    

if nPar ~= length(cases)
    error('Nr of cases <<%d>> must equal number of ParArrays <<%d>>',...
          length(cases),nPar);
end

Pars = cell(nVar,nPar);

if isempty(cases)
    for ip=1:size(Pars,2)-2
        cases{ip} = sprintf('case %d',ip);
    end
end

for iv=1:nVar
    varName = varargin{1}{iv,1};
    Pars{iv,1} = varName;              % col1: variabe name
    Pars{iv,2} = varargin{1}{iv,end};  % col2: use of this variable

    for ip = 1:nPar
        % first assert variable names are equal
        if ~strcmp(varName,varargin{ip}{iv,1})
            error('varname %s ~= varname in input %d',varName,ip);
        end
        % then add value to cell array
        Pars{iv,ip+2} = varargin{ip}{iv,2};
    end
end

        
fprintf('%12s%5s','Parameter','use');
fprintfs('%12s',cases);
fprintf('%12s','stdPar');
fprintf('%12s','%uncert');
fprintf('\n');

j=0;
for iv = 1:nVar
    fprintf('%12s%5d',Pars{iv,1},Pars{iv,2});
    for ip = 3:size(Pars,2)
        fprintf('%12.4g',Pars{iv,ip});
    end
    if Pars{iv,2}
        j=j+1;
        fprintf('%12.4g',stdP(j));
        fprintf('%12.4g',uncertP(j));
    end
    fprintf('\n');
end

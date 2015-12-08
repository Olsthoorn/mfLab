function X = arrayCheck(name,X)
%ARRAYCHECK verfies that there are no NaN's or Infs in the array.
%
% USAGE:
%   X = arrayCheck(arrayName,X)
%
%   Verfies that there are no NaN's or Infs in the array.
%   Used in mf_setup to prevent layer crashes of MODFLOW and related
%   programs when excuting. Causes of such crashes may be hard to find if
%   they are due to sporadic NaN in one of the input files.
%
%   The reason is, that if there, are Modflow doesn't see this and it
%   may take a long time to figure out what the problem is.
%
%   TO 130301

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

    if iscell(X) 
        for i=1:numel(X)
            check(name,X{i});
        end
    else
        check(name,X);
    end

end

function check(name,X)
    if any(isnan(X(:)))
        error('%s: Your array %s has one or more NaNs. USGS programs can''t handle that.',mfilename,name);
    end

    if any(X(:)==Inf)
            error('%s: Your array %s has one or more Inf. USGS programs can''t handle that.',mfilename,name);
    end

    if any(X(:)==-Inf)
            error('%s: Your array %s has one or more -Inf. USGS programs can''t handle that.',mfilename,name);
    end
end
function done(toc)
%DONE write done with of without the toc value
%
% USAGE:
%    done(toc)
%
% Convenient to insert in scripts or functions to show that something that
% takes time has finished and how long it took
%
% TO 120615

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin==0
    fprintf('done\n');
else
    frintf('done, %.2f seconds\n',toc);
end
function B = mf_setTime(basename,B)
%MF_SETTIME adds simulation time to budget struct
%
% USAGE:
%    B = mf_setTime(B)
% 

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

sp_time = sp_timeObj(basename);

for it = 1:numel(B)
    B(it).time = sp_time.t(B(it).period,B(it).tstp);
end
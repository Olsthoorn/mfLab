function ok = BCN_Data_Check(type)
%BCN_DATA_CHECK check if workspace has data for this BCN (Boundary Condition Info)
%
% USAGE:
%    ok = BCN_Data_Check(type,who) -- check if workspace has data for this BCN
%
%    type is one of ['WEL', 'GHB', 'RIV', 'DRN', 'CHD' ...]
%    who is the cell array with all variables names in the workspace
%
% Example:
%    WEL = [1 3 2 1 24.0];
%    ok=BCN_Data_Check('WEL');
%
%
% TO 120528

ok=strmatchi(type,who);

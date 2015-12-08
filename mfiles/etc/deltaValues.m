function C = deltaValues(C,STCONC)
%DELTAVALUES substract starting values from data in struct which must have field 'values'
%
% USAGE:
%    C = deltaConc(C,STCONC)
%
% used in animateObj and animate functions
%
% TO 130403

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

for i = 1:numel(C)
    C(i).values = C(i).values-STCONC;
end
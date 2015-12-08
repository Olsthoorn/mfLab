function y=ierfc(z,n)
%IERFC repeated integral of complementary error function
%
% USAGE:
%     y=ierfc(z,n)
%
% needed in transient 1D groundwater analysis. I.e. in the Edelman and
% Bruggeman (1999) series.
%
% See Stegun and Abramowitz (1964) Handbook of mathematical functions,
% Item 7.2.7, p299
%
% TO 100222

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin<2
    error('%s: not enough input arguments',mfilename);
end

n=round(n);

if n<-1 
    error('second argument n must be integer >= -1');
end
switch n
    case -1
        y=2/sqrt(pi)*exp(-z.^2); return;
    case 0
        y=erfc(z); return;
    otherwise
        y=-z/n.*ierfc(z,n-1)+(1/2/n)*ierfc(z,n-2);
end

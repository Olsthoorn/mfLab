function range=binRange(x,dx)
%BINRANGE gets a range min(x):dx:max(x) for the data x
%
% USAGE:
%    range=binRange(x,dx)
%
% dx is chosen range interval. The number of significant didgets is also obtained
% from dx
%
% TO 100529 100811

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

m=min(x(:));
M=max(x(:));

N=(M-m)/dx+2;

if m<0, n=fix(m/dx)-1; else n=fix(m/dx); end

range=dx*(n:n+N);

if length(range)>250
    
    error(['binRange length must be < 200, check input\n',...
          'Lowest=%g step=%g Highest=%g, number of contours: %d\n'],m,dx,M,length(range));
end

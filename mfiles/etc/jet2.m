function cmap=jet2(len,I)
%JET2 generate useful colormap similar to jet
%
% USAGE:
%    cmap=jet(len,I)
%
% gives color I if I is specified (?)
%
% generates colormap used for salinity figures
%
% TO 110320 110416

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin<2, len=64; end
    
cmap =[linspace(0,1,len)',...
      [linspace(0,1,len/2)'; linspace(1,0,len/2)'],...
       linspace(1,0,len)'...
      ];

  colormap(cmap);
  
  if exist('I','var');
      cmap=cmap(I,:);
  end

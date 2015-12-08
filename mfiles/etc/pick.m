function picked=pick(it,values,dim)
%PICK picks next value from list and start with 1  after reaching end
%
% USAGE:
%    picked=pick(it,values[,dim]); 
%
% Examples:
%     clr=pick(it,'brgkwyb');
%     val=pick(it,[3 5 12 4 35]);
%
% TO 120101 120531 121126

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin<2, picked=[]; return; end

if nargin<3
    NC = numel(values);
else
    NC = size(values,dim);
end

ic=rem(it,NC); ic(ic==0)=NC;

if iscell(values)
    picked = values{ic};
elseif nargin<3
    picked = values(ic);
else
    switch dim
        case 1
            picked = values(ic,:,:);
        case 2
            picked = values(:,ic,:);
        case 3
            picked = values(:,:,ic);
        otherwise
            error('%s: dim can be 1,2 or 3',mfilename);
    end
end

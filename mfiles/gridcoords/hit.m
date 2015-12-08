function Idx = hit(xGr,a,b)
%HIT get index of cells between a and b or if one input argument the cell in which a resides
%
% Example:
%   Idx = hitting(xGr,a,b)
%
% See also: between after before inside outside beyond above below inMesh
%
%   TO 1204020

if exist('xGr','var') && isempty(xGr),
    Idx=[];
    return;
else
    xGr = squeeze(xGr);
end

switch nargin
    case {0 1}, Idx=NaN; return;
    case 2, b=a(1); a=a(1);
    otherwise 
       if b<a(1), dum=b(1); b=a(1); a=dum; end
end
   

   yGr=1:length(xGr);
   
   I = floor(interp1(xGr(:),yGr,[a b]));
   
   Idx=min(I):max(I);
   Idx=Idx(~isnan(Idx));
   
   
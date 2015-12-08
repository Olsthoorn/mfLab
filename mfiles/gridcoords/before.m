function Idx = before(xGr,a)
%BEFORE gets indices of xGr or xm before point a, where xGr is ascending
%
% Example:
%    Idx = before(xGr,a)
%
% See also:  hit between fallsIn above below after beyond outside inMesh
%
%   TO 1204020

   yGr=1:length(xGr);
   ym =0.5*(yGr(1:end-1)+yGr(2:end));
      
   ix = interp1(xGr(:),yGr,a);
   
   Idx=find(ym<ix);
   
   if isempty(Idx)
       Nx=length(ym);
       if xGr(end)>xGr(1)
           if a>xGr(end)
               Idx=1:Nx;
           end
       else
           if a<xGr(end)
               Idx=1:Nx;
           end
       end
   end
      
     
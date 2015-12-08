function [conc,totim]=plotCNC(co1,co2,CNC,period,tstp,varargin) % idx,C,clr,contfun)
%PLOTCNC plots the concentration contained in struct CNC for given stress
%   period and time step on plane given by dir and index idx
%   co1 and co2 are the cell center coordinates (xm,ym or xm,zm or ym,zm)
%
%   dir is either 'x','y','z' or 'row','col','layer', indiating the plane to contour
%   idx is the corresponding plane number
%   C are countours, clr is the colors and contfun the handle to a contour
%   function @contour or @contourf
%   colorEge may be used to color the edges in contourf, set to 'nono' to
%   remove the edges. Leave out for default or use one of 'none' 'r' 'b' 'g' 'k' etc
%
%EXMPLE
%   [conc,totim]=plotCNC(co1,co2,CNC,period,tstp[,dir[,idx[,C[,clr[,contfun[,colorEdge]]]]]])
%
% TO 090222

% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later


L=length(varargin);

if L<1 || isempty(varargin{1}), dir='row';  else dir=varargin{1};  end
if L<2 || isempty(varargin{2}), idx=1;      else idx=varargin{2};  end
if L<3 || isempty(varargin{3}), C='';       else C  =varargin{3};  end
if L<4 || isempty(varargin{4}), clr='';     else clr=varargin{4};  end
if L<5 || isempty(varargin{5}), contfun=@contourf; else contfun=varargin{5}; end
if L<6 || isempty(varargin{6}), colorEdge=NaN; else colorEdge=varargin{6}; end

[conc totim]=plotHDS(co1,co2,CNC,period,tstp,dir,idx,C,clr,contfun,colorEdge);

function [hds totim]=plotHDS(co1,co2,H,period,tstp,varargin) % idx,C,clr,contfun)
% [hds totim]=plotHDS(co1,co2,H,period,tstp[,dir[,idx[,C[,clr[,contfun[,colorEdge]]]]]])
% plots the heads or concentrations contained in struct H for given
% stress period and time step on plane given by dir and index idx
% co1 and co2 are the cell center coordinates (xm,ym or xm,zm or ym,zm)
% dir is either 'x','y','z' or 'row','col','layer', indiating the plane to contour
% idx is the corresponding plane number
% C are countours, clr is the colors and contfun the handle to a contour
% function @contour or @contourf
% colorEdge is optional, may be used in contourf to remove edgecolors
% between colored regions (leave out for default or use one of 'none' 'r' 'b' 'g' 'k' etc
% hds are the heads as they are plotted, i.e. permuted
% totim is the total simulation time of these heads
% TO 090222


% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later


L=length(varargin);

if L<1 || isempty(varargin{1}), dir='row';  else dir=varargin{1};  end
if L<2 || isempty(varargin{2}), idx=1;      else idx=varargin{2};  end
if L<3 || isempty(varargin{3}), C=NaN;      else C  =varargin{3};  end
if L<4 || isempty(varargin{4}), clr=NaN;    else clr=varargin{4};  end
if L<5 || isempty(varargin{5}), contfun=@contourf; else contfun=varargin{5}; end
if L<6 || isempty(varargin{6}), colorEdge=NaN; else colorEdge=varargin{6}; end

dir=lower(dir(1));  % only use lowercase of first letter of dir

for i=1:length(H)
    found=H(i).period==period && H(i).tstp==tstp;
    if found,
        if isfield(H,'totim'), totim=H(i).totim; end
        if isfield(H,'time'),  totim=H(i).time;  end
        break;
    end
end
if ~found, error('Can''t find period %d and tstp %d in struct',period,tstp); end

switch dir
    case {'y' 'r'}, hds=permute(H(i).values(idx,:,:),[3 2 1]);
    case {'x' 'c'}, hds=permute(H(i).values(:,idx,:),[3,1,2]);
    case {'z' 'l'}, hds=permute(H(i).values(:,:,idx),[1 2 3]);
    otherwise
        error('direction %c illegal, used ''row'' ''col'' or ''layer''',dir);
end

if isnan(C)
    [c,hdl]=contfun(co1,co2,hds);
elseif isnan(clr)
    [c,hdl]=contfun(co1,co2,hds,C);
else
    [c,hdl]=contfun(co1,co2,hds,C,clr);
end
if ~isnan(colorEdge), set(hdl,'edgecolor',colorEdge); end
hold on

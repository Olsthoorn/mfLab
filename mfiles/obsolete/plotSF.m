function sf=plotSF(co1,co2,B,period,tstp,dir,idx,psi,clr,contfun)
% Psi=plotSF(x,z,B,period,tstp[,dir[,idx[,psi[,clr[,contfun]]]]])
% plots the stream function contained in budget struct B for given
% stress period and time step on x-section given by dir and index idx
% co1 and co2 are the grid coordinates
% dir is one of 'x','y','z','row','col','layer' indicating plane to contour
% x and row are equivalent as are y and col as are z and layer
% default is row
% idx is the corresponding plane number, default = 1
% psi are countours, clr are colors default is empty
% contfun is the handle to a contour functio contfun=@contour or @contourf or any
% default is contour
% legal alternative

if ~exist('dir','var') || isempty(dir), dir='row';    end
if ~exist('idx','var') || isempty(idx), idx=1;        end
if ~exist('psi','var') || isempty(psi), psi=NaN;      end
if ~exist('clr','var') || isempty(clr), clr=NaN;      end
if ~exist('contfun','var') || isempty(contfun), contfun=contour; end

dir=lower(dir(1));  % only use lowercase of first letter of dir

for i=1:length(B)
    found=(B(i).period==period && B(i).tstp==tstp);
    if found, break, end
end
if ~found,     error('Can''t find timestep %d of stress period %d',tstp,period); end

switch dir
    case {'y' 'r'}, LBL='FLOWR'; % plot xz plane (row) for given iRow
    case {'x' 'c'}, LBL='FLOWF'; % plot yz plane (column) for given iCol
    case {'z' 'l'}, LBL='FLOWL'; % plot xy plane (layer) for given iLay
    otherwise
        error('use COL ROW or LAYER as direction');
end

for jLBL=1:length(B(i).label)
    found=findstr(LBL,B(i).label{jLBL});
    if found, break; end
end
if ~found,  error('Can''t find %s in budget file, check direction in call',LBL); end
        
switch dir
    case {'y' 'r'}, sf=permute(B(i).term{jLBL}(idx,:,:),[3 2 1]);
    case {'x' 'c'}, sf=permute(B(i).term{jLBL}(:,idx,:),[3,1,2]);
    case {'z' 'l'}, sf=permute(B(i).term{jLBL}(:,:,idx),[1 2 3]);
end
sf=[flipud(cumsum(flipud(sf)));...
    zeros(1,size(sf,2))];  % stream function

if isnan(psi)
    psi=50;  % just 50 psi contour lines
elseif length(psi)==1  % treat psi as dPsi
        psi=sort([-psi:-psi:min(sf(:)), 0, psi:psi:max(sf(:))]);
end

if isnan(psi)
    contfun(co1(2:end-1),co2,sf(:,1:end-1));
elseif isnan(clr)
    contfun(co1(2:end-1),co2,sf(:,1:end-1),psi);
else
    contfun(co1(2:end-1),co2,sf(:,1:end-1),psi,'color',clr);
end
hold on

function total = showTotal(o,varargin)
%GRIDOBJ/SHOWTOTAL: -- contours the total of gr.DZ.*(sum arg(1),arg(2), ...),3)
% using about 50 contour lines.
%
% USAGE: kDtotal = gr.showTotal([ax,]HK[n,contourOptions]);
% HK is a full array
% n is number of contours, defalt is 50
% contourOptions are options accepted by contour (matlab function)
%
% the usage example computes and plots tot total transmissivity of the model
%
% TO 131011

k=0;

[ax,varargin] = getType(varargin,'axis',[]);
if isempty(ax)
  k=1;
  figure('name','total see title axis','pos',screenPos(0.75));
  ax = axes('nextplot','add');
  xlabel(ax,'x [m]');
  ylabel(ax,'y [m]');
end

total = ones(o.size);
nIn = numel(varargin);
s = '';
for i=1:nIn
    [A,varargin] = getNext(varargin,'double',[]);
    if isempty(A)
        break
    end
    if numel(A)==1
        continue
    end
    total = total .* A;
    k=k+1;
    s = [s '+' inputname(k)]; %#ok
end

s=s(2:end);

total = sum(total .* o.DZ,3);

title([' sum ( gr.DZ .* sum(' s ',3),3)']);

range = contRange(total);

% compute decent contour line range
% M = max(total(:)); m=min(total(:)); dT = (M-m)/nContour;
% d = 10^floor(log10(dT)); if dT/d >5, d=5*d; elseif dT/d>2, d=2*d; end
% M=d*(floor(M/d)+1); m = d*(floor(m/d));
% range = m:d:M;

[C,h] = o.contourf(ax,total,range); clabel(C,h,'color','w');

hb=colorbar; set(get(hb,'title'),'string',['sum(gr.DZ.*(' s '),3)']);

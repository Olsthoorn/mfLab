function h = blockInterp(h,xgr,xm)
% interpolate h onto xm if h is stepwise constant in cells mplied by xGr
% neatly handle NaNs in h
% TO 120803

DELTAX = 0.1;

hds = [h(1) h(1:end-1);h(2:end) h(end)];

% remove NaN's from fist row of h=[hds,h] only for symmetrical results
% if nargin<4
%     
%     if isnan(hds(1)), hds(1)=find(~isnan(hds'),1); end % make sure left head ~= NaN, arbitrary is ok
% 
%     I=find(isnan(hds));
%     while ~isempty(I)
%         hds(I)=hds(I-1); % take left neighbor
%         I=find(isnan(hds));
%     end
% end
%
% hds = [hds; h];

xgr = [xgr(1:end-1)-DELTAX; xgr(2:end)+DELTAX];

h   = interp1(xgr(~isnan(hds)),hds(~isnan(hds)),xm);

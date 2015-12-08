function Rxx = crosscorr(x,y,lags)
% Rxx = crosscorr(x,y,lags) -- computes cross corr. between vectors x and y
% x and y must be of same size.
% TO 141009

if nargin<3, lags = 25; end

x=x(:); N=numel(x);

if nargin<2
    y=x;
else
    y=y(:);
    if ~all([numel(x),numel(y)]==N)
        error('inputs x and y must be of same size');
    end
end

Rxx = zeros(1,lags);
for i=1:lags
    C = corrcoef(x,y);
    Rxx(i) = C(end,1);
    y = [y(2:end); y(1)];
end

function Rxx = autocorr(x,lags)
% Rxx= aucorr(x,lags) :-- computes autocorrelation of vector x
% SEE ALSO crosscorr
% TO 141009

if nargin<2, lags=25; end

Rxx = crosscorr(x,x,lags);


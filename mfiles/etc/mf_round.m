function x = mf_round(x,N)
%MF_ROUND round to n decimals (only needed for Matlab version < R2014)
% For later versions just use round(x,n)
% TO 151112

eN = 10.^N;
x = round(x * eN) / eN;
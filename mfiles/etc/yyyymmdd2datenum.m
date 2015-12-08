function D = yyyymmdd2datenum(D)
% Transform YYYYMMDD to datenum (YYYYMMDD is used in KNMI data files)
%
% TO 141009

YR  = floor(D/10000);
MON = round(100*(floor(D/100)/100 - YR));
DAY = round(100*(D/100 - floor(D/100)));
D   = datenum(YR,MON,DAY);

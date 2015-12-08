function KNMIdata()
%% Help for getting and importing KNMI data files
%
% KNMI is the Royal Netherlands Meteorological Institute in De Bilt,
% Netherlands.
% Among all it does, it provides data of meteostations measured throuhout
% the past. These data are available on their site
%

web('http://www.knmi.nl/klimatologie/monv/reeksen/')

% for daily precipitation data of over 300 stations
% and on

web('http://www.knmi.nl/klimatologie/daggegevens/download.html');

% for daily data of all measured parameters for about 30 stations.
%
% This site also gives a map of the stations, their number and their name.
%
% Hourly and daily data of all elements are available from 1901 for 30
% stations.
% Precipitation data, daily values, are available since 1906 for 300
% stations.
% These are probably the most useful data for time series analysis.
%
% ET values are generally available on the 30 stations providing all
% measured items, but only since 1987. This is the so-called Makkink ET, i.e.
% the reference ET computed using only temperature and short wave incoming
% radiation. It is valid for well watered grass land. It is generally used
% to compute the ET for well watered crops in the Netherlands throug
% multiplication by a crop factor that depends both on the crop type and on
% the month in the year. These crop factors can be found in the KNMI info
% file
%
% monv_toelichting.pdf
%
% which is also in this folder (in Dutch).
%
% The data files can be downloade as zip files that yield text files with
% comma separated data after unzipping. Header lines provide information on
% the contents of the columns. Values unmeasured are blank.
%
% The text files of the 30 stations with all parameters have name
%
% etmgeg_???.txt
%
% where ??? is the station number.
%
% The text files of the 300 stations with only precipiation have name
%
% neerslaggeg_!!!!!_????.txt
%
% where !!!! is the name of the station and ??? its
% number.
%
% mfLab has to mfiles to import KNMI data directly from the text files
%
% TPE = KNMIimport_etmgeg(stationNr)
% TP  = KNMIimport_neerslaggeg(stationNr)
%
% both accept a station number as input
% TPE is an array with three columns [datenum, precip [m/d], Emak [m//d]]
% TP  is an array with two   columns [datenum, precip [m/d] ]
%
% notice that empty values for P and E in the data files become NaN
%
% Also notice that the datenum is the time at which the precipitation was
% registered, which is 08:00 UT (universal time, i.e. 9:00h in winter and
% 10:00h in summer. The taime covered is from 08:00 on the previous day to 08:00
% on the current day (see file monv_toelichting.pdf). 
% The ET values are for the day registered from 00:00 to 23:59:59

help KNMIdata

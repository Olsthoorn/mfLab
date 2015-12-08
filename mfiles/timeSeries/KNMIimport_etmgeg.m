function [TPE,stationName] = KNMIimport_etmgeg(varargin)
%% Read KNMI data file -- read KNMI data txt file
%
% Read ASCII file for general Meteo stations of KNMI
%
% USAGE: 
%    [TPE,stationName] = KNMIimport_etmgeg(station [,KNMIdataPath] [,'NaNok'])
%    [TPE,stationName] = KNMIimport_etmgeg(station [,'.'[,'NaNok'])  % current directory
%    [TPE,stationName] = KNMIimport_etmgeg(station [,'NaNok']) % data on default directory
%
% station       = KNMI station number
% KNMIdataPath  = path to directory with KNMIdata files
% TPE           = [datenum precip[m/d] Emakkink[m/d]]
% stationName   = name of the meteo station, obtained from list
%                 meteoStations.csv in default KNMI directory, see next
%
% default KNMI data directory:
%       'mflab/examples/TimeSeriesAnalysis/KNMI-data-files/'
%
% The data originate from
%
% 'http://www.knmi.nl/klimatologie/daggegevens/download.html'
%
% The data is one of the 30 KNMI stations with all data in a file name
%
% 'etmgege_???.txt'
%
% where ??? is the KNMI meteo station number. E.g. 240 is Schiphol.
%
% See also: KNMIimport_neerslaggeg  KNMIdata
%
% TO141021

% Hard wired path with in the mflab directory structure

    [NaNok,varargin] = getWord(varargin,'NaN');

    [station,varargin] = getNext(varargin,'numeric',[]);
    if nargout>1 || ~isempty(station)
        [Path,~] = getNext(varargin,'char',[]);
        if isempty(Path)
            Path = 'mflab/examples/TimeSeriesAnalysis/KNMI-data-files/';
            pwdPath = pwd;
            i       = strfind(pwdPath,'mflab');
            if ~i
                error(['your current path is not in the mflab directory structure\n',...
                    'Use path to KNMI data files as second argument of call.']);
            end
            Path = [pwdPath(1:i-1) Path];
        end
        Path = fileparts(Path);
        KNMIdataPath = Path;
    end
        % get number of header lines by scanning for # in first column
    if ~isempty(station)
        try
            fname = sprintf([P 'etmgeg_%d.txt'],station);
        catch %#ok
            fname = sprintf('./etmgeg_%d.txt' ,station);
        end
    else
        fname = getNext(varargin,'char',[]);
        d = dir(fname);
        if numel(d)>1
            error('file name must be unique %s',fname);
        elseif numel(d)<1
            error('can''t find file %s',fname);
        end
    end

    fid = fopen(fname,'r');
    if fid<0    
        web('http://www.knmi.nl/klimatologie/daggegevens/download.html');
        error(['Can''t find file, get your file at the site shown.' fname]);
    end

    nHeaderLines = lookfor(fid,'# STN',0);

    fclose(fid);    

    A         = importdata(fname,',',nHeaderLines);
    hdr       = strtrim(A.textdata(end,:));
    I = strmatchi({'YYYYMMDD','RH','EV24'},hdr,'exact');

    T = A.data(:,I(1)); PE = A.data(:,I(2:3));

    PE(PE(:,1)<0,1) = 0;  % values -1 mean <0.05 mm
    PE=PE/1e4;            % from 0.1 mm/d to m/d

    t = datenum(floor(T(:,1)/10000),...
                floor(T(:,1)/100)-100*floor(T(:,1)/10000),...
                T(:,1)-100*floor(T(:,1)/100)) + 8/24;

    TPE = [t PE];

    if ~NaNok
        TPE = TPE(~ (isnan(TPE(:,2)) | isnan(TPE(:,3))),:);
    end

    if nargout>1
        station = regexp(fname,'_([0-9]+)','tokens');
        station = sscanf(station{1}{1},'%f');

        %Path = fileparts(Path);

        % get meteo station name from KNMI default diretory
        stations    = importdata(fullfile(KNMIdataPath,'meteoStations.csv'),',',1);
        i           = find(stations.data(:,1) == station);
        stationName = stations.textdata{i+1,1};
    end

    function n = lookfor(fid,word,n1)
         for i=1:1000 %#ok
             s = fgetl(fid);
             %fprintf('%s\n',s)
             if numel(s)>numel(word) && strcmp(word,s(1:numel(word)));
                 n = n1+i;
                 return;
             end
         end
         error('word <<%s>> not foud in this file',word);
    end
end


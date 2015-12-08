function [TP,stationName] = KNMIimport_neerslaggeg(varargin)
%% Read KNMI data file -- read KNMI data txt file
%
% Read ASCII file for preciptiation Meteo stations of KNMI (there are 300
% such stations in the Netherlands)
%
% USAGE: 
%    [TP,stationName] = KNMIimport_neerslaggeg(fname [,'NaNok'])
%    [TP,stationName] = KNMIimport_neerslaggeg(fname [,'NaNok])  % default directory
%    [TP,stationName] = KNMIimport_neerslaggeg(station[,path['NaNok]]);
%
% where fname may contain full path plus wild cards but it should pertain
% to a single, unique data file fname
%
% station       = KNMI station number
% stationName   = KNMI station name from file neerslagStations in KNMI
% default directory (see below)
%
% KNMIdataPath  = path to directory with KNMIdata files
% TP            = [datenum precip[m/d]]
%
% default KNMI data directory:
%       'mflab/examples/TimeSeriesAnalysis/KNMI-data-files/'
%
% The data originate from
%
% 'http://www.knmi.nl/klimatologie/daggegevens/download.html'
%
% The data is one of the 300 KNMI stations with only precipitation having
% in a file name
%
% 'neerslaggeg_!!!!!!_???.txt'
%
% where !!!!! is the name of the stationand
% ??? is the KNMI meteo station number. E.g. 069 is Appelscha.
%
% See also: KNMIimport_etmgeg  KNMIdata
%
% TO141021

    [NaNok,varargin] = getWord(varargin,'NaN');

    [station,varargin] = getNext(varargin,'numeric',[]);
    
    if nargout>1|| ~isempty(station)

        [Path,~] = getNext(varargin,'char');
        if isempty(Path)
            % Default KNMI data path, will be generated from current
            % location
            Path = 'mflab/examples/TimeSeriesAnalysis/KNMI-data-files/';
            pwdPath = pwd; % Current path, assumed to be on the mflab directory structure

            i = strfind(pwdPath,'mflab');
            if ~i
                error(['your current path is not in the mflab directory structure\n',...
                    'Use path to KNMI data files as second argument of call.']);
            end
            Path = [pwdPath(1:i-1) Path];
        end
    end

    if ~isempty(station)
        
        % does the file with this station number exist?
        try
            d = dir([Path sprintf('neerslaggeg*%d.txt',station)]);
            fname = [Path     d.name];
        catch %#ok try current directory
            d = dir(['./' sprintf('neerslaggeg*%d.txt',station)]);
            fname =  ['./' d.name];
        end
    else    
        [fname,~] = getNext(varargin,'char',[]);
        d = dir(fname);
        if numel(d)>1
            error('specified file must be unique');
        elseif numel(d)<1
            error('can''t find file %s\n',fname);
        end
    end

    fid = fopen(fname,'r');
    if fid<0
        web('http://www.knmi.nl/klimatologie/monv/reeksen/')
        error(['Can''t find file %s. Get your file from this shown website' fname]);
    end

    nHeaderLines = lookfor(fid,'STN',0);
    nHeaderLines = lookfor(fid,'STN',nHeaderLines);

    fclose(fid);   

    A         = importdata(fname,',',nHeaderLines);
    I = strmatchi({'YYYYMMDD','RD'},strtrim(A.colheaders),'exact');

    T = A.data(:,I(1)); P = A.data(:,I(2));

    P(P(:,1)<0,1) = 0;  % values -1 mean <0.05 mm
    P=P/1e4;            % from 0.1 mm/d to m/d


    t = datenum(floor(T(:,1)/10000),...
                floor(T(:,1)/100)-100*floor(T(:,1)/10000),...
                T(:,1)-100*floor(T(:,1)/100)) + 8/24;

    TP = [t P];

    if ~NaNok
        TP = TP(~isnan(TP(:,2)),:);
    end

    % get station name from KNMI default diretory
    if nargout>1
        station = regexp(fname,'_([0-9]+)','tokens');
        station = sscanf(station{1}{1},'%f');

        Path = fileparts(Path);
        stations    = importdata(fullfile(Path,'neerslagStations.csv'),',',1);
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

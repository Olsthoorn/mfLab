classdef meteoObj
%METEOOBJ --- data from KNMI historic database
% BRON:
% KONINKLIJK NEDERLANDS METEOROLOGISCH INSTITUUT (KNMI)
% 
% SOURCE:
% ROYAL NETHERLANDS METEOROLOGICAL INSTITUTE (KNMI)
% 
% TO 130519
    properties        
        hdr
        data
        descr
        date
    end
    methods
        function o=meteoObj(meteoFile)
            %METEOOBJ -- constructor meteoObj
            %
            % USAGE:
            %   mObj = meteoOjb(meteoFileName)
            % INPUT
            %   meteoFileName, from historic database of www.KNMI.nl
            %   for a single meteo station.
            %
            % Example
            %     see usage
            %
            % TO 130519
            
            if nargin<1
                return;
            end
            
            fid = fopen(meteoFile,'r'); if fid<0, error('%s: Can''t open file <<%s>>',mfilename,meteoFile); end

            while true
                sHdr=fgets(fid);
                if sHdr(1)=='#'
                    break;
                end
            end
            
            sHdr(sHdr=='#')=' ';
            sHdr(sHdr==',')=' ';
            o.hdr={};
            while true
                [o.hdr{end+1},sHdr] = strtok(sHdr); %#ok
                if isempty(sHdr)
                    break;
                end
            end
            
            if isempty(o.hdr{end}), o.hdr = o.hdr(1:end-1); end

            fgets(fid);
                   
            f0     = ftell(fid);
            s      = fgets(fid);
            f1     = ftell(fid);
            reclen = f1-f0;
            fseek(fid,0,1);
            fe     = ftell(fid);
            Nrow   = (fe-f0)/reclen;
            Ncol   = numel(strfind(s,','));
            format = repmat(' %d',[1,Ncol]);
            fseek(fid,f0,-1);

            o.data = textscan(fid,format,Nrow,'delimiter',',','collectOutput',1);
            o.data = o.data{1};

            fclose(fid);
            
            % get first EV24
            EVT=o.data(:,strmatchi('EV24',o.hdr));
            i = find(EVT>0,1,'first');
            o.data = o.data(i:end,:);
            
            o = o.setDescr();
            o = o.setDate();
        end
        function o = setDescr(o)
            %DESCR -- set description of headers
            o.descr.YYYYMMDD  =  'Datum (YYYY=jaar MM=maand DD=dag) / Date (YYYY=year MM=month DD=day)';
            o.descr.DDVEC     =  'Vectorgemiddelde windrichting in graden (360=noord, 90=oost, 180=zuid, 270=west, 0=windstil/variabel). Zie http://www.knmi.nl/klimatologie/achtergrondinformatie/windroos.pdf / Vector mean wind direction in degrees (360=north, 90=east, 180=south, 270=west, 0=calm/variable)';
            o.descr.FHVEC     =  'Vectorgemiddelde windsnelheid (in 0.1 m/s). Zie http://www.knmi.nl/klimatologie/achtergrondinformatie/beaufortschaal.pdf / Vector mean windspeed (in 0.1 m/s)';
            o.descr.FG        =  'Etmaalgemiddelde windsnelheid (in 0.1 m/s) / Daily mean windspeed (in 0.1 m/s) ';
            o.descr.FHX       =  'Hoogste uurgemiddelde windsnelheid (in 0.1 m/s) / Maximum hourly mean windspeed (in 0.1 m/s)';
            o.descr.FHXH      =  'Uurvak waarin FHX is gemeten / Hourly division in which FHX was measured';
            o.descr.FHN       =  'Laagste uurgemiddelde windsnelheid (in 0.1 m/s) / Minimum hourly mean windspeed (in 0.1 m/s)';
            o.descr.FHNH      =  'Uurvak waarin FHN is gemeten / Hourly division in which FHN was measured';
            o.descr.FXX       =  'Hoogste windstoot (in 0.1 m/s) / Maximum wind gust (in 0.1 m/s)';
            o.descr.FXXH      =  'Uurvak waarin FXX is gemeten / Hourly division in which FXX was measured';
            o.descr.TG        =  'Etmaalgemiddelde temperatuur (in 0.1 graden Celsius) / Daily mean temperature in (0.1 degrees Celsius)';
            o.descr.TN        =  'Minimum temperatuur (in 0.1 graden Celsius) / Minimum temperature (in 0.1 degrees Celsius)';
            o.descr.TNH       =  'Uurvak waarin TN is gemeten / Hourly division in which TN was measured';
            o.descr.TX        =  'Maximum temperatuur (in 0.1 graden Celsius) / Maximum temperature (in 0.1 degrees Celsius)';
            o.descr.TXH       =  'Uurvak waarin TX is gemeten / Hourly division in which TX was measured';
            o.descr.T10N      =  'Minimum temperatuur op 10 cm hoogte (in 0.1 graden Celsius) / Minimum temperature at 10 cm above surface (in 0.1 degrees Celsius)';
            o.descr.T10NH     =  '6-uurs tijdvak waarin T10N is gemeten / 6-hourly division in which T10N was measured; 6=0-6 UT, 12=6-12 UT, 18=12-18 UT, 24=18-24 UT ';
            o.descr.SQ        =  'Zonneschijnduur (in 0.1 uur) berekend uit de globale straling (-1 voor <0.05 uur) / Sunshine duration (in 0.1 hour) calculated from global radiation (-1 for <0.05 hour)';
            o.descr.SP        =  'Percentage van de langst mogelijke zonneschijnduur / Percentage of maximum potential sunshine duration';
            o.descr.Q         =  'Globale straling (in J/cm2) / Global radiation (in J/cm2)';
            o.descr.DR        =  'Duur van de neerslag (in 0.1 uur) / Precipitation duration (in 0.1 hour)';
            o.descr.RH        =  'Etmaalsom van de neerslag (in 0.1 mm) (-1 voor <0.05 mm) / Daily precipitation amount (in 0.1 mm) (-1 for <0.05 mm)';
            o.descr.RHX       =  'Hoogste uursom van de neerslag (in 0.1 mm) (-1 voor <0.05 mm) / Maximum hourly precipitation amount (in 0.1 mm) (-1 for <0.05 mm)';
            o.descr.RHXH      =  'Uurvak waarin RHX is gemeten / Hourly division in which RHX was measured';
            o.descr.PG        =  'Etmaalgemiddelde luchtdruk herleid tot zeeniveau (in 0.1 hPa) berekend uit 24 uurwaarden / Daily mean sea level pressure (in 0.1 hPa) calculated from 24 hourly values';
            o.descr.PX        =  'Hoogste uurwaarde van de luchtdruk herleid tot zeeniveau (in 0.1 hPa) / Maximum hourly sea level pressure (in 0.1 hPa)';
            o.descr.PXH       =  'Uurvak waarin PX is gemeten / Hourly division in which PX was measured';
            o.descr.PN        =  'Laagste uurwaarde van de luchtdruk herleid tot zeeniveau (in 0.1 hPa) / Minimum hourly sea level pressure (in 0.1 hPa)';
            o.descr.PNH       =  'Uurvak waarin PN is gemeten / Hourly division in which PN was measured';
            o.descr.VVN       =  'Minimum opgetreden zicht / Minimum visibility; 0: <100 m, 1:100-200 m, 2:200-300 m,..., 49:4900-5000 m, 50:5-6 km, 56:6-7 km, 57:7-8 km,..., 79:29-30 km, 80:30-35 km, 81:35-40 km,..., 89: >70 km)';
            o.descr.VVNH      =  'Uurvak waarin VVN is gemeten / Hourly division in which VVN was measured';
            o.descr.VVX       =  'Maximum opgetreden zicht / Maximum visibility; 0: <100 m, 1:100-200 m, 2:200-300 m,..., 49:4900-5000 m, 50:5-6 km, 56:6-7 km, 57:7-8 km,..., 79:29-30 km, 80:30-35 km, 81:35-40 km,..., 89: >70 km)';
            o.descr.VVXH      =  'Uurvak waarin VVX is gemeten / Hourly division in which VVX was measured';
            o.descr.NG        =  'Etmaalgemiddelde bewolking (bedekkingsgraad van de bovenlucht in achtsten, 9=bovenlucht onzichtbaar) / Mean daily cloud cover (in octants, 9=sky invisible)';
            o.descr.UG        =  'Etmaalgemiddelde relatieve vochtigheid (in procenten) / Daily mean relative atmospheric humidity (in percents)';
            o.descr.UX        =  'Maximale relatieve vochtigheid (in procenten) / Maximum relative atmospheric humidity (in percents)';
            o.descr.UXH       =  'Uurvak waarin UX is gemeten / Hourly division in which UX was measured';
            o.descr.UN        =  'Minimale relatieve vochtigheid (in procenten) / Minimum relative atmospheric humidity (in percents)';
            o.descr.UNH       =  'Uurvak waarin UN is gemeten / Hourly division in which UN was measured';
            o.descr.EV24      =  'Referentiegewasverdamping (Makkink) (in 0.1 mm) /  Potential evapotranspiration (Makkink) (in 0.1 mm)';
        end
        function plot(o,varargin)
            %PLOT -- plots daily precipitaiton and evapotranspiration as
            %afunction of time
            [ax     ,varargin] = getProp(varargin,'axis',[]);
            [pos    ,varargin] = getProp(varargin,{'figPos','pos'},[]);
            [figName,   ~    ] = getProp(varargin,'fig',[]);
            
            if ~isempty(figName)
                if isempty(pos), pos = screenPos(0.75); end
                figure('name',figName,'position',pos);
                ax = axes('nextplot','add');
            else
                if isempty(ax)
                    ax=gca;
                    set(ax,'nextplot','add');
                end
            end
            xlabel('time'); ylabel('mm');
            RH   = o.data(:,strmatchi('RH'  ,o.hdr,'exact'));
            EV24 = o.data(:,strmatchi('EV24',o.hdr,'exact'));
            bar(ax,o.date,[RH -EV24]/10); % bacause in 10th of mm
            datetick(ax);
            %bar(o.date,EV24,'r');
        end
        function o = setDate(o)
            %SETDATE -- convert time vector into datenum
            ddVec = double(o.data(:,strmatchi('YYYYMMDD',o.hdr)));
            Y = floor(ddVec/10000);
            M = floor(ddVec/100);
            M = 100 * (M/100-floor(M/100));
            D = 100*(ddVec/100-floor(ddVec/100));
            o.date = datenum(round(Y),round(M),round(D));
        end
    end
end
    
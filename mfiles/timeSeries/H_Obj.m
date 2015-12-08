classdef H_Obj < IN_Obj
    % Menthanthes H_Obj, holding an explaining time series like
    % PREC or EVAP. The H_Obj is a child of IN_Obj
    % Under construction TO 141027
    properties
        tnocode
        measpointlev
        date
        handmeas
        vegtype
        area
        areacode
        diver_files
        comment
    end
    methods
        function o = H_Obj(varargin)
            % H_Obj -- constructor of menynathes H_Obj (Head files)
            %
            % USAGE: headobjects = H_Obj([P *.csv]);
            %
            % where P is the path to the directory with the tno head data
            % csv files
            %
            % headobjects = array of H_Obj holding the data of each tno
            % head data file
            %
            % TO 141027
            
            o = o@IN_Obj(varargin{:});
            % get the data of head series
            
            if nargin==0, return; end
            
            [P,~,~] = fileparts(varargin{1});
            
            d = dir(varargin{1});
            Nfiles = numel(d);
            fprintf('reading H files ...\n');
            for io=Nfiles:-1:1
                F = fullfile(P,d(io).name);
                O(io) = o.readTNOgwst(F);
                
                 % assume measured at noon
                O(io).values(:,1) = floor(O(io).values(:,1))+0.5;
                               
                fprintf('.'); if rem(Nfiles-io+1,50)==0, fprintf('%d\n',Nfiles-io+1); end
            end 
            if rem(Nfiles-io+1,50)~=0, fprintf('%d',Nfiles-io+1); end
            fprintf(' H-files read.\n');
            o = O;
        end
        function o = readTNOgwst(o,fname)
            % H_Obj/readTNOgwst(fname) ---read a TNO head data file
            % These file are downloaded from Dino-loket.nl and stored as
            % csv text files
            % 
            % USAGE o = readTNOgwst(fname);
            %
            % TO 141027

            fp =  fopen(fname,'r');

            %% Info in TNO groundwater head files
            % One record per line:
            % Titel:,,,,,,,,,,,
            % Gebruikersnaam:,,,,,,,,,,,
            % Periode aangevraagd:,01-01-1800,tot:,11-09-2014,,,,,,,,
            % Gegevens beschikbaar:,27-01-1961,tot:,31-12-2010,,,,,,,,
            % Datum: ,11-09-2014,,,,,,,,,,
            % Referentie:,NAP,,,,,,,,,,
            % 
            % NAP:,Normaal Amsterdams Peil,,,,,,,,,,
            % MV:,Maaiveld,,,,,,,,,,
            % MP:,Meetpunt,,,,,,,,,,
            % 
            % Locatie,Filternummer,Externe aanduiding,X-coordinaat,Y-coordinaat,Maaiveld (cm t.o.v. NAP),Datum maaiveld gemeten,Startdatum,Einddatum,Meetpunt (cm t.o.v. NAP),Meetpunt (cm t.o.v. MV),Bovenkant filter (cm t.o.v. NAP),Onderkant filter (cm t.o.v. NAP)
            % B06D0110,001,06DP0110,190573,578183,70,13-07-1971,25-07-1996,14-01-2003,29,-41,-1293,-1393
            % B06D0110,001,06DP0110,190573,578183,55,14-01-2003,14-01-2003,31-12-2010,49,-6,-1273,-1373
            % 

            %% Info rgarding elevation and its changes over time
            % Locatie
            % Filternummber
            % Externe_aanduiding
            % xcoord      =
            % ycoord      =
            % mvNAP       =
            % mvdatum     = 
            % startdatum  =
            % einddatem   =
            % mpntNAP     =
            % mpntMV      =
            % bkfNAP      =
            % okfNAP      =

            o.comment    = fname;
            o.type       = 'head';

              date_      = lookfor(fp,'Datum');
            o.date       = ['TNO file date= ' date_{end}];

            Hdr          = lookfor(fp,'Locatie'); %#ok

            % Read the data that specifies elevation and the period over
            % which the elevation is valid (one line per period)
            for iL=1:100
                s = fgetl(fp);
                %fprintf('%s\n',s);
                if isempty(s), break; end
                C = dlmread(s,',');
                o.UserData.mvDatum(iL,1)      = datenum(flipud(abs(C{7}))');
                o.UserData.startDate(iL,1)    = datenum(flipud(abs(C{8}))');
                o.UserData.endDate(iL,1)      = datenum(flipud(abs(C{9}))');
                o.UserData.mpntNAP(iL,1)      = C{10}/100;
                o.UserData.mpntMV(iL,1)       = C{11}/100;
                o.UserData.bkfNAP(iL,1)       = C{12}/100;
                o.UserData.okfNAP(iL,1)       = C{13}/100;
            end
            %% Store data in menyanthes struct
            o.name       = sprintf('%s_%d',C{ 1},C{2});
            o.tnocode    = C{ 1};
            o.filtnr     = C{ 2};
            o.xcoord     = C{ 4};
            o.ycoord     = C{ 5};
            o.surflev    = o.UserData.mpntNAP(end);  % in m not in cm
            o.upfiltlev  = o.UserData.bkfNAP(end);
            o.lowfiltlev = o.UserData.okfNAP(end);

            %% Non-Menyanthes info --> UserData
            o.UserData.extName      = C{3};

            %% Actual data of the time series in the file

            % Get header:
            % Locatie,Filternummer,Peildatum,Stand (cm t.o.v. MP),Stand (cm t.o.v. MV),Stand (cm t.o.v. NAP),Bijzonderheid,Opmerking,,,
            Hdr = lookfor(fp,'Locatie');
            o.UserData.DataHdr = Hdr(2:end); % discard repeated name of tube

            % Actual data: 
            % Locatie,Filternummer,Peildatum,Stand (cm t.o.v. MP),Stand (cm t.o.v. MV),Stand (cm t.o.v. NAP),Bijzonderheid,Opmerking,,,
            % B06D0110,001,27-01-1961,118,,,,,,,,
            % B06D0110,001,27-02-1961,114,,,,,,,,
                        
            % Read the data from the remainder of the tno file
            C = textscan(fp,'%s %f %f-%f-%f %f %f %f %s %s %f %f %f %f %f',Inf,'delimiter',',','emptyValue',NaN);

            % Get the numeric data preceded by date:
            o.values = [datenum(C{5},C{4},C{3}),  [C{6}, C{7}, C{8}]/100 ];

            % Also get remarks
            o.UserData.remark  = {o.values(:,1) C{9} C{10}};

            fclose(fp);

            function tokens = lookfor(fp,word)
                % tokens -- split the header lines properly into tokens
                % but first look for the line starting with a given word
                for ii=1:100
                    s = fgetl(fp);  %fprintf('%s\n',s);
                    n = length(word);
                    if length(s)>=n && strcmpi(s(1:n),word)
                        [tokens] = regexp(s, '[A-Za-z0-9\_\+\-\.\(\) ]*', 'match');         
                        return;        
                    end
                end
                fprintf('h''m can''t find token\n');
            end
            function token = dlmread(s,delimiter)
                I=[0 find(s==delimiter) length(s)+1];
                for j=numel(I)-1:-1:1;
                    txt = s(I(j)+1:I(j+1)-1);
                    token{j} = sscanf(txt,'%s');
                    if isempty(token{j})
                        token{j}=NaN;
                    else
                        token{j} = sscanf(txt,'%f');
                        if isempty(token{j})
                            token{j}=txt;
                        end

                    end
                end
            end
        end
        function plot(o,varargin)
            % H_Obj/plot -- plot the H_Obj array
            %
            % USAGE: H.plot([lineSpec])
            %
            % TO 141028
            
            [ax,varargin] = getType(varargin,'axis',gca);
            set(ax,'nextPlot','add','xgrid','on','ygrid','on');
            if isempty(varargin)
                for io=1:numel(o)
                    plot(ax,o(io).values(:,1),o(io).values(:,end),mf_color(io));
                end                
            else
                for io=1:numel(o)
                    plot(ax,o(io).values(:,1),o(io).values(:,end),varargin{:});
                end
            end
            legend({o.name});
            datetick();
        end
        function mdl = calibrate(H,fun,pIn,IN,dtSampling)
            % Calibrate -- calibrate the time series
            %
            % USAGE: H = H.calibrate(fun,pIn,IN);
            %
            % TO 141028
            
            if nargin<4
                dtSampling = round(mean(diff(IN(1).values(:,1))));
            end
            
            mdl = tsaMdlObj(fun,pIn,H,IN);
            
            [t,P,E] = IN.prepareData(dtSampling);

            mdl = mdl.calibrate(t,P,E); % note o(io)=H
        end
        function prepareData(o)
            error('prepareData only works on IN_Obj class objects not on class %s',class(o));
        end
    end
end
    
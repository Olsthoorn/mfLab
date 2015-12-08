classdef IN_Obj % -- Menyanthes explaining time series (Extension .IN)
    properties
        name       = 'Name of the explaining series';
        type       = 'none';  % 'PREC','EVAP'
        values     = NaN;
        xcoord     = NaN;
        ycoord     = NaN;
        surflev    = NaN;
        filtnr     =  NaN;
        upfiltlev  = NaN;
        lowfiltlev = NaN;
        datlog_serial = NaN;
        R          = NaN;
        UserData
    end
    methods
        function o=IN_Obj(type,varargin)
            % Constructor of Menyanthes explaining time series (IN_Obj)
            % read data from a KNMI time series, which can be a KNMI meteo
            % station with all measured variables or a KNMI precipitation
            % station with only daily precipitation data.            
            %
            % USAGE: IN_Obj = IN_Obj(type,KNMIstationNr)
            % USAGE: IN_Obj = IN_Obj(type,KNMIfilesdirectory))
            %
            % type = 'PREC' | 'EVAP
            %
            % See also: KNMIdata KNMIimport_etmgeg, KNMIimport_neerslaggeg
            %
            % For IN obj see Menyanthes files or manual
            %
            % TO 141027
            
            if nargin==0, return; end
            
            if ~strcmpi(class(o),'IN_Obj'), return; end
            
            if nargin <2
                error('USAGE: IN_obj = IN_Obj(type,KNMMfiles)');
            end
                        
            [P,~,~] = fileparts(varargin{1});
            
            d = dir(varargin{1}); Nfiles=numel(d);
            fprintf('reading IN files ...\n');
            for io=Nfiles:-1:1
                O(io) = o.KNMIimport(type,fullfile(P,d(io).name));
                fprintf('.');
                if rem(Nfiles-io+1,50)==0, fprintf('%d\n',Nfiles-io+1); end
            end 
            if rem(Nfiles-io+1,50)~=0, fprintf('%d',Nfiles-io+1); end
            fprintf(' done !\n');
            o = O;
        end
        function plot(o,varargin)
            % IN_Obj/plot --- plots the data of the IN_Obj
            %
            % USAGE:  IN_Obj.plot(varargin)
            % where varargin are options as valid for plot in Matlab
            %
            % TO 141027
            
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
        function [R,o] = distance(o,x,y)
            % Obj/distance -- compute distance R between point and well|gauge
            %
            % USAGE: [R,H_Obj] = H_Obj.distance(x,y);
            %
            % TO 141028
            
            R = sqrt((x-[o.xcoord]).^2+(y-[o.ycoord]).^2);
            if nargout>1
                for io=numel(o):-1:1
                    o(io).R = R(io);
                end
            end
        end
        function plotLoc(o,varargin)
            % Obj/plotLoc -- plot location of gauges or wells
            %
            % USAGE: Obj.plotLoc([plotCodes]);
            %
            % TO 141028
            
            [ax,varargin] = getType(varargin,'axis',gca);
            set(ax,'nextPlot','add','xGrid','on','yGrid','on');
            
            if isempty(varargin)
                plot(ax,[o.xcoord],[o.ycoord],'o');
            else
                plot(ax,[o.xcoord],[o.ycoord],varargin{:});
            end
        end
        function o = KNMIimport(o,type,fname)
            o.type = type;
            switch o.type
                case 'PREC'
                    try
                        [o.values,o.name] = KNMIimport_neerslaggeg(fname,'NaNok');
                         o.values         = o.values(:,1:2);
                    catch ME
                        [o.values,o.name] = KNMIimport_etmgeg(fname,'NaNok');
                        fprintf('%s\n',ME.message);
                    end
                case 'EVAP'
                    [o.values,o.name] = KNMIimport_etmgeg(fname,'NaNok');
                case 'WEL'
                    error('type %s not yet implemented',o.type);
                case 'RIV'
                    error('type %s not yet implemented',o.type);
                otherwise
                    error('unknown type %s: use PREC or EVAP',o.type);
            end
        end
        function o = cleanup(o,fac)
            % cleanUp (remove NaN in data (heads in NAP))
            % if fac given, remove outliers, i.e. mean+/-fac*std
            %
            % USAGE: Obj.cleanup([fac]);
            %
            % TO 141104
            for io=numel(o):-1:1
                o(io).values = o(io).values(~isnan(o(io).values(:,end)),:);
                if nargin>1 && fac>0
                    mu = mean(o(io).values(:,end));
                    s  = std( o(io).values(:,end));
                    J  = o(io).values(:,end) < mu + fac*s & ...
                         o(io).values(:,end) > mu - fac*s;
                     o(io).values = o(io).values(J,:);
                end
            end
        end
        function o = period(o,T,N)
            % H_Obj/selectPeriod -- select wells with data in period
            %
            % USAGE: period = Obj.period();  % shows time span of data
            %           Obj = Obj.period([t1 t2]); truncate to period [t1 t2]
            %           Obj = Obj.period([t1 t2],N); same but at least N
            %        measurements within period
            %
            % TO 141028
            
            if nargin==1
                for io=1:numel(o)
                    fprintf('%s  %s\n',...
                        datestr(o(io).values(  1,1)),...
                        datestr(o(io).values(end,1)));
                end
                return;
            else
                for io=numel(o):-1:1
                    t = o(io).values;
                    o(io).values = o(io).values(t>=min(T) & t<=max(T),:);                
                end
            
                % Then a minimum of N data should be present
                if nargin>2
                    for io=numel(o):-1:1
                        if size(o(io).values,1)<N
                            o(io)=[];
                        end
                    end
                end
            end
        end
        function [tSampling,P,E] = prepareData(o,dtSampling)
            % Prepare IN data, used by simulate
            % Works only on IN_Obj
            
            P = o(strmatchi('PREC',{o.type}));
            E = o(strmatchi('EVAP',{o.type}));

            if ~isempty(P)
                tp =  P(1).values(:,1);
                P  =  P(1).values(:,end);
                tp = tp(~isnan(P));
                P  =  P(~isnan(P));
            else
                tp = [];
            end
            
            if ~isempty(E)
                te =  E(1).values(:,  1);
                E  =  E(1).values(:,end);
                te = te(~isnan(E));
                E  =  E(~isnan(E));
            else
                te = [];
            end
            
            if  isempty(tp)
                if sempty(te)
                    error('no input data series with type PREC and or EVAP are given');
                else
                    tp = te;
                    P = zeros(size(tp));
                end
            else
                if isempty(te)
                    te = tp;
                    E = zeros(size(te));
                else
                    % skip
                end
            end
            
            % match times
            tMin = max(tp(  1),te(  1));
            tMax = min(tp(end),te(end));
            P    = P( tp>=tMin & tp<=tMax);
            tp   = tp(tp>=tMin & tp<=tMax);
            E    = E( te>=tMin & te<=tMax);
            te   = te(te>=tMin & te<=tMax);
            
            if nargin<2, dtSampling=1; end
            
            tSampling = (tp(1):dtSampling:tp(end))';
            % this is not correct
            
            dtp = diff(tp);
            dte = diff(te);
            
            warning('off','all'); % NaN will be resolved separately
            P         = interp1(tp,cumsum(P.*[dtp(1); dtp]),tSampling);
            E         = interp1(te,cumsum(E.*[dte(1); dte]),tSampling);
            P         = [0; diff(P)/dtSampling];
            E         = [0; diff(E)/dtSampling];
            warning('on','all');
        end
    end
end
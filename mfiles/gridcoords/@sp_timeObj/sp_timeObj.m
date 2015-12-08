classdef sp_timeObj
    %SP_TIMEOBJ class def of timeObj to get time for every stress period and time step
    %
    %  sp_timeObj holds the time at the end of every stress period and
    %  timestep of a MODFLOW or MT3MDS/SEAWAT simulation.
    %  fields have the names of come columns of the PER worksheet used in
    %  the workbooks associated with mfLab models.
    %
    % USAGE:
    %    o = sp_timeObj(basename[,tzero]);
    %    one object will be generated for each stress period.
    %    the number of output objects equals the number of stress periods.
    %
    %    tzero can be a datenum allowing to relat times to actual dates
    %
    % What it does:
    %    gets a sp_timeObj containing the time at the end of every
    %    time step of every stress period. It is computed from the
    %    data in the PER worksheet of the workbook associated with
    %    each mfLab example.
    %    This sp_timeObj is useful to compute the times with every
    %    element of the Budget, which only caries time step and stress
    %    period, but not the actual time.
    %    Methods will be provied to extract specific information.

    % TO 130401
    
    properties
        iPer  % [scalar] stress period nr
        nStp  % [salar] nr of time steps in stress period obj
        iStp  % [vector] nrs of time steps within stress period [1...NStp(SP(i))]

        tsMult % [scalar] time step multiplier in this stress period obj

        perlen   % [scalar] length of this stress period in days (assumed to be time unit of model)
        spStart  % [scalar] absolute starting time (datenum) of this stress period
        spEnd    % [scalar] absolute ending   time (datenum) of this stress period

        tsLen    % [1,NStp] length of time steps in this stress period in days (assumed to be time unit of model)
        tsStart  % [1,NStp] absolute starting time (datenum) of time steps in this stress period  
        tsEnd    % [1,NStp] absolute ending   time (datenum) of time steps in this stress period

        t0 = 0;  % [scalar] absolute time of start of simulation (may be set to zero)
    end
    properties (Dependent=true)
        spLen   % same as perlen
    end
    methods
        function o = sp_timeObj(XLSF,varargin)
        %SP_TIMEOBJ constructor of sp_timeObj class
        %
        % USAGE:
        %    o = sp_timeObj(basename[,tzero]);
        %    one object will be generated for each stress period.
        %    the number of output objects equals the number of stress periods.
        %
        %    tzero can be a datenum allowing to relat times to actual dates
        %
        % What it does:
        %    gets a sp_timeObj containing the time at the end of every
        %    time step of every stress period. It is computed from the
        %    data in the PER worksheet of the workbook associated with
        %    each mfLab example.
        %    This sp_timeObj is useful to compute the times with every
        %    element of the Budget, which only caries time step and stress
        %    period, but not the actual time.
        %    Methods will be provied to extract specific information.
        %
        % TO 130404
            
            if nargin<1, return; end

            [tzero, ~ ] = getNext(varargin,'double',0);

            % required headers from the PER sheet
            Hdrs = {'iPER','PERLEN','NSTP','TSMULT'};
            
            % get these headers and their values, plus number of stress
            % periods
            [hdr,PERdata,NPER] = getPeriods(XLSF,'PER',Hdrs);

            % determine column numbers associated with certain data
            IPER   = PERdata(:,strmatchi('iPER'  ,hdr));
            PERLEN = PERdata(:,strmatchi('PERLEN',hdr));
            NSTP   = PERdata(:,strmatchi('NSTP'  ,hdr));
            TSMULT = PERdata(:,strmatchi('TSMULT',hdr));

           [o(1:NPER).t0] = deal(tzero);
           
           for io=1:NPER
               o(io).iPer       = PERdata(IPER(io));               
               o(io).perlen     = PERLEN(io);
               o(io).nStp       = NSTP(io);
               o(io).tsMult     = TSMULT(io);
           end
           
           previous = tzero;
           for io=1:NPER
                o(io).spStart = previous;
                o(io).spEnd   = previous+o(io).perlen;

                % power series based on TSMULT
                cp = cumprod(ones(1,o(io).nStp)*o(io).tsMult);

                o(io).tsLen   = o(io).perlen * cp/sum(cp);

                ts =  [0 cumsum(o(io).tsLen)];
                o(io).iStp    = 1:o(io).nStp;
                o(io).tsStart = previous + ts(1:end-1);
                o(io).tsEnd   = previous + ts(2:end);

                previous    = o(io).spEnd;
           end
        end

        function spLen = get.perlen(o)
        %SP_TIMEOBJ/PERLEN -- gets stress period length
            spLen = o.perlen;
        end
        
        function o = sett0(o,tzero)
        %SP_TIMEOBJ/T0 --- sets absolute datenum
            for io=1:numel(o)
                o(io).spStart = o(io).spStart - o(io).t0 + tzero;
                o(io).spEnd   = o(io).spEnd   - o(io).t0 + tzero;
                o(io).tsStart = o(io).tsStart - o(io).t0 + tzero;
                o(io).tsEnd   = o(io).tsEnd   - o(io).t0 + tzero;
                o(io).t0      = tzero;
            end
        end

        function [Hdr,Lst] = list(o)
            %SP_TIMEOBJ/LIST --- generates the list [SP STEP time]
            Hdr = {'SP', 'TSTP', 'TIME'};
            N = sum([o.nStp]);
            Lst = NaN(N,3);
            k=1;
            for i=1:numel(o)
                for j=1:o(i).nStp
                    Lst(k,:) = [o(i).iPer, o(i).iStp(j), o(i).tsEnd(j)];
                    k=k+1;
                end
            end
            % remove what's empty
            Lst(k:end,:)=[];                
        end
        
        function show(o)
            %SP_TIMEOBJ/SHOW --- displays the object
            [Hdr,Lst] = o.list();
            fprintf('%7s  %7s  %7s',Hdr{:}); fprintf('\n');
            fprintf('%5d  %5d %12g\n',Lst');
        end
        
        function tme = t(o,SP,TSP)
            %SP_TIMEOBJ/T --- yields datenum for end of given SP and time step TSP
            % time = timeObj.t;          --> times at end of stress periods
            % time = timeObj.t(n)        --> time at end of stress period SP
            % time = timeObj.t(o,SP,TSP) --> time at end of tstep TSP in stress period SP            
            switch nargin
                case 1
                    tme = [o.spEnd] + o(1).t0;
                case 2
                    tme = o(SP).tsEnd(end)+o(SP).t0;
                otherwise
                    tme = o(SP).tsEnd(TSP)+o(SP).t0;
            end
        end
    end
end
classdef wellSeriesObj
    properties
        nr; % well series nr or id (see id) depending on use
        id; % for instance a GIS_id or other pointer to a well series
        wellType;
        name;
        longName;
        wellCodes;
        nrOfWells;
        layer;
        shdl=[];  % handel to well series graphics object        
        remark    % for any purpose
        UserData; % for any purpose
        created   % time stamp   
        Dt; t;    % stress period length and time at end of stres period
        Q; C; Cout; NCOMP; %Plugged in in mf_setWellSeries
        children; % ids or nrs of the wells belonging to this wells series        
    end
    properties (Dependent = true)
        numberOfWells
    end
    methods
        function o=wellSeriesObj(well,parentLabel)
            % wellSeries = wellSeriesObj(well,parentName)
            % generates well series objects.
            % parentName is the field in well.UserData which holds the
            % wellSeries group number to which the well belongs. Make sure
            % you have that field in the well specification table from
            % which the wellObj are geneated.
            % TO 130307

            if nargin==0,
                return;
 
            end

           wellNrs = [well.nr    ]; % the well ids in the order of the wells in the series
           
           if nargin>1               
               if ~ischar(parentLabel)
                   error('Second argument must be name of a field in object or in its UserData');
               end
               fnms = fieldnames(well);
               if strmatchi(parentLabel,fnms)
                   for iw = numel(well):-1:1
                       well(iw).parent = well(iw).(parentLabel);
                   end
               elseif ~isempty(well(1).UserData)
                   fnms= fieldnames(well(1).UserData);
                   if strmatchi(parentLabel,fnms)
                       for iw=numel(well):-1:1
                           well(iw).parent = well(iw).UserData.(parentLabel);
                       end
                   else
                       error('well.UserData has not field <<%s>>',parentLabel);
                   end
               else
                   eror('well has no field <<%s>>',parentLabel);
               end
           end
           
           parents = [well.parent];

           
            if strcmp(class(well),'MNW1Obj') || strcmp(class(well),'wellObj')
                %well is here an array of MNW
                wellSeriesId = unique([well.parent]);
                
                for is = numel(wellSeriesId):-1:1
                    o(is).wellType = class(well);
                    o(is).nr       = wellSeriesId(is);
                    o(is).id       = o(is).nr;
                    kids = find(parents == o(is).nr);
                    o(is).children = wellNrs(kids); 
                    o(is).name =  well(kids(1)).name;
                    for ic = numel(o(is).children)-1:-1:1
                        if ~strcmpi(o(is).name,well(kids(ic)).name)
                            o(is).name =  sprintf('%s%d','wellSeries',o(is).nr);
                            break;
                        end
                    end
                    
                    %% Get Children for this well series
                    IW  = find([well.parent] == o(is).nr);
                    for i=numel(IW):-1:1
                        o(is).children(i) = well(IW(i)).nr;
                    end
                    o(is).created  = now;
                    o(is).shdl = []; % handle to well series graphics obj
                end
                fprintf('ready\n');
            else
                error('%s: first argument of call must be of clas <<wellObj>> or <<MNW1Obj>> not <<%s>>',mfilename,class(well));
            end
        end
        function numberOfWells = get.numberOfWells(o), numberOfWells = numel(o.children); end        
        
        function stamp  =stamp(o,fmt), stamp=datestr(o.created,fmt); end
        function created=get.created(o), created=o.created; end
                
        function o=setCout(o,C,iComp), o=setWellCout(o,C,iComp); end
        function Cout=get.Cout(o), if isempty(o.Cout), Cout=NaN(size(o.C)); else Cout=o.Cout; end; end
        
        
        function plotQ(   o      ,varargin), plot(o.t,o.Q  ,varargin{:}); end
        function plotC(   o,iComp,varargin), plot(o.t,iComp,varargin{:}); end
        
        function plotCout(o,iComp,varargin), plot(o.t,o.Cout(iComp,:),varargin{:}); end
        function plotCin( o,iComp,varargin), wellPlotC(o.t,iComp,varargin{:}); end
        
        function plotXY(o,well,varargin)
            % wellSeries.plotXY(well,varargin);
            % plots location of wells in wellSeries on map
            % 130307
            wellNrs = [well.nr];
            for iw = find(ismember(o(iws).children,wellNrs))
              well(iw).plotXY(varargin{:});
            end
        end
        function plot3D(o,well,varargin)
            % wellSeries.plot3D(well,varargin)
            % 3D plot of the wells in the well serie 
            % 130307
            wellNrs = [well.nr];
            for iw = find(ismember(o(iws).children,wellNrs))
                well(iw).plot3D(varargin{:});
            end
        end
        function well=plot2D(o,well,varargin)
            % wellseriesObj.plot2D(well,iPer)  refresh well status for this stress period.
            % wellseriesObj.plot2D(well,options) as usual: 'facecolor'.'y', 'edgecolor','k'
            % wellseriesObj.plot2D(well,'y',options) same but uses well's ycoordinate instead of x
            wellNrs = [well.nr];
            for iw = find(ismember(o(iws).children,wellNrs))
                well(iw) = well(iw).plot2D(varargin{:});
            end
        end
    end
end

classdef observationObj
    %OBSERVATIONBOBJ class definition of observation objects
    %
    % USAGE:
    %    see explantion at observationObj/observatinObj (constructor below)
    %
    % observationObj defines observations for use in groundwater models. The objects hold
    % property values and share methods that facilitate using observations.
    % Some methods need a grid object. But a grid object is not necessary
    % to create well objects. This is because upon creation, well objects
    % only contain their pertinent data. Later on, when a grid object is
    % available they may be put into the grid, upon which additional grid
    % values are added to the observations, so that each well then also holds its
    % position in the grid. See the methods of observationObj and of gridObj for
    % more information.
    %
    % SEE ALSO gridObj wellObj
    %
    % TO 120807 130620
    properties
        nr; id;      % well nr or id (see id)
        parent;   % well series id to which well belongs (wellSeries(..).id)
        name=''; longname='';
        legend=''; % legend text, if empty, name is used

        x; y; z;
        ix; iy; iLay; idx;
        LRC       % Layer row col for the screen penetrated cells
        DZ;       % thickness of mode layers penentrated by screen
        T;        % transmissivity of layer
        fQ;       % fraction of total Q from each model layer
        Dt; t;    % stress period length and time at end of stres period
        NCOMP=0;  % C=specified injection conc,
        species   % names of the species/components
        lineSpec  % lineSpec for plotting (default mf_color(iw))
        lineWidth=1;
        wpix=2;   % pixelwidth of observations on screen for plotting        
        whdl=NaN(3,1);  % 3-valued handle: well casing and well screen in XSdrawing
                  %  whdl(1)=screen
                  %  whdl(2)=casing
                  %  whdl(3)=marker on casing top        
        
        remark = ''; % any text
        code,     % for whatever purpose
        created   % moment at which well was create by the programm
        FaceColor = 'w';  % outline color of plotted observations in XS or YS
        EdgeColor = 'b';  % screen coor of plotted observations in XS or YS        
 
        UserData  % for any purpose, contains the actual observation data
                  % see o(iw).UserData.H(it).value etc.
                  %     o(iw).UserData.DDN(it).value
                  %     o(iw).UserData.(speciesName).value
    end
    methods
        function o=observationObj(varargin)
            %OBSERVATIONOBJ constructor of observation objects
            % observation objects are almost the same as wells as they are
            % in fact wells that do not have an extraction. Their calling
            % syntax is therefore also compapitible. In the future, they
            % may be totally the same, execpt for their class name.
            %
            % USAGE:
            %   obs = observationObj(basename,sheetNm);
            %   obs = observationObj(basename,sheetNm,gr[,HK]);
            %   obs = observationObj(basename,sheetNm,gr[,HK],{species,quantities});
            %   
            % First call, generate the observation objects from the data
            % provided in the sheetNm of the workbook basename. The only
            % requirements are to provide nr,x,y,z1,z2.
            %
            % To put them into an existing grid, supply a gridObj, here
            % indicated with gr.
            % To allow the concentration to be computed as a weighted
            % average based on the transmissivities of the penetrated
            % model layers also supply the full 3D HK array (horizontal k).
            % To load data into the observation supply the names of the
            % species in the case of concentrations (or temperature). The
            % names are immaterial, but the given names will be attributed
            % to the first and second etc. conc file MT3D00n.UCN (where n
            % is a sequence number (1,2,3 ...). To load head, drawdown etc
            % specify the specific works 'head' and or 'drawdown' in the
            % quantaties input. Place the species and quantities in a cell
            % array as {'chloride','temp','head'}. The observation will
            % then read the data for the first and second concenration in
            % the simulation and calls them 'chloride' and 'temp' in that
            % order. It further regonizes the word 'head' and will load the
            % file [basename '.HDS'] and extract the heads at the
            % obsservation point location from it.
            %         'grid',gridObj
            %% To make it work:
            %   put a list of nr,x,y,z1,z2 points in sheet sheetNm of workbook
            %   basename.xls. These points determine the observation point
            %   locations in 3D.
            %   Then call the constructor.
            % 
            % EXAMPLE:
            %   piezoms = observationObj(basename,sheetNm,gr,,HK,'head');
            %   piezoms.plot('head')          % plots head of all piezoms
            %   piezomes([1 4]).plot('temp')  % plots only for nrs.1 and 4
            %
            %   gr must be of class gridObj. It may be placed anywhere in
            %   the input or be given as parameter, value pair 'grid',gr
            %
            %   As is the case with wellObj, HK is desired to compute a
            %   transmissivity weighted average of the obervation wells
            %   along their screen if it passes more than one model layer.
            %
            % EXAMPLE:
            %   obsWells = observationObj(basename,sheetNm,gr[,HK],{'salinity','temp','head'});
            %   obsWells = observationObj(basename,sheetNm,gr[,HK],{'salinity','temp','drawdown'});
            %   obsWells = observationObj(basename,sheetNm,gr[,HK],{'salinity','temp','heads'},...
            %       'STCONC',STCONC,'STRTHD',STRTHD);
            %   obsWells.plot('salinity',varargin)
            %   obsWells.plot('temp',varargin);
            %   obsWells.plot('head',varargin);
            %   obsWells.plot('drawdown',varargin)
            %
            %   varargin may be used to supply property,value pair
            %   information for the plot (same options as with ordinary
            %   plot function).
            %
            %  'head' is recognized, it will load the heads from [basename'.HDS']
            %  'drawdown' is recognized and will be loaded from [baseanem '.DDN']
            %   This file must exist if it is used this way.
            %
            %   Cell array {'salinity','temp','head'} will also be
            %   recognized. The word 'head' is always interpreted as heads and
            %   'drawdown' as drawdown.
            %
            %   The other strings are interpreted as the names of the two
            %   species in the simulation. NCOMP must then be 2.
            %   The concentration for these two species will be read from
            %   files MT3D00n.UNC with n the species number.
            %
            %   The property value pair 'STRTHD',STRTHD may be put in the input
            %   call. In that case the drawdown will be calculated using
            %   [basename '.HDS'] and STRTHD.
            %
            %   The property value pair 'STCONC',STCONC may be put in the input
            %   call. In that case the change of concentration will be
            %   computed by subtracting the STCONC for the species from its
            %   simulated values.
            %
            % SEE ALSO: gridObj.well (shortcut)
            %
            % TO 120807
            
            if nargin<1,
                return;
            end            
            
            %% grid? See if you can get a grid object from the input 'grid',gr
            [gr,varargin] = getProp(varargin,'grid',[]);            
      
            [nxyz,varargin] = getNext(varargin,'double',[]);
            if ~isempty(nxyz)
                if size(nxyz,2)~=4
                    error('%s: list if obs points must be of shape [n x y z]',mfilename);
                else
                    USAGE = 1;
                end
            else                
                [basename,varargin] = getNext(varargin,'char',[]);
                [sheetNm ,varargin] = getNext(varargin,'char',[]);
                if ~isempty(basename) && ~isempty(sheetNm)
                    USAGE = 2;
                else
                    error('%s: illegal ussage see <<help %s>>',mfilename,mfilename);
                end
            end

            if USAGE == 1
                for i=size(nxyz,1):-1:1
                    o(i).nr= nxyz(i,1);
                    o(i).id= nxyz(i,1);
                    o(i).name = sprintf('Obs%03d',nxyz(i,1));
                    o(i).longname = o(i).name; o(i).legend = o(i).name;
                    o(i).x = nxyz(i,2);
                    o(i).y = nxyz(i,3);
                    o(i).z = nxyz(i,4);
                    o(i).created=now;
                    o(i).lineSpec = mf_color(i);
                end
                o = o([o.nr]>0);
                return;
            end

            % Continue with USAGE 2
            members = fieldnames(o);
            table   = getTable(basename,sheetNm,members,'Hor');
            field   = fieldnames(table);

            for j=1:numel(field)
                for iw=numel(table):-1:1
                    o(iw).lineSpec   = mf_color(iw); % default lineSpec
                    o(iw).(field{j}) = table(iw).(field{j});
                end
            end

            % if legend (legend text) is not specified, use name
            for iw=numel(o):-1:1
                if isempty(o(iw).legend)
                    o(iw).legend = o(iw).name;
                end
            end
            
            %% deal with possible fields like z1 and z2 separately,
            % zScreen..., that all join into screen z coordinates.
            % TODO
            %   This is meant for curved, i.e. multi section screens and may need some work and
            %   precision to be robust
            if isfield(table,'UserData')
                    UDfield = fieldnames(table(1).UserData);  % UDfield = User Data field names
            end
            
            if isempty(o(iw).z)
                L = regexpis(UDfield,'z[1-9]|zscr'); % indices of fields z1 z2 z3 z4 ... etc
                if isempty(L)
                    error('%s: Add columns z1 z2 to well sepcification table to specify top and bottom of well screen',mfilename);
                else
                    for iw=numel(o):-1:1
                        for i=length(L):-1:1
                            o(iw).z(i)=o(iw).UserData.(UDfield{L(i)});
                        end
                    end
                end
            end

            % clean up if necessary
            for iw=numel(o):-1:1
                o(iw).x(isnan(o(iw).x))=[];
                o(iw).y(isnan(o(iw).y))=[];
                o(iw).z(isnan(o(iw).z))=[];
            end


            % if not see if grid is given as one of the inputs
            if isempty(gr) && ~isempty(varargin)
                i = strmatchi('gridObj',cellfun(@class,varargin,'uniformOutput',false));
                if i(1)
                    gr = varargin{i(1)};
                    varargin(i)=[];
                else                    
                    gr = [];
                end
            end
            
            % put observation into grid
            if isempty(gr)
                return;
            else
                [HK,varargin] = getNext(varargin,'double',[]);

                if isempty(HK)
                    o = o.toGrid(gr);
                else
                    % To deal with grid consisting of a single layer for which
                    % matlab drops the 3rd dimension in size we use grSize and
                    % arSize in the check for compabible input
                    grSize= gr.size;
                    arSize= size(HK);

                    if ~isnumeric(HK) || ~all(arSize==grSize(1:numel(arSize)))
                        error(['%s: Fourth argument of call must be a numerical array of size <<%d %d %d>>\n',...
                            'You should provide the model array of TRAN or HK\n',...
                            'You provided and object of class <<%s>>\b',...
                            'REMEDY: Check if grid has LAYCBD compatible iwth HK or TRAN'],...
                            mfilename,gr.size,class(varargin{1}));
                    end
                    o = o.toGrid(gr,HK);  % toGrid(gr,HK) (original HK)
                end
            end

            %% At this point all observations have their idx, so we can continu
            
            %% STRTHD requested?
            [STRTHD,varargin] = getProp(varargin,'STRTHD',[]);
            if ~isempty(STRTHD)
                for io = numel(o):-1:1
                    o(io).UserData.STRTHD = STRTHD(o(io).idx);
                end
            end
            %% STCONC requested?
            [STCONC,varargin] = getProp(varargin,'STCONC',[]);
            if ~isempty(STCONC)
                if ~iscell(STCONC), STCONC={STCONC}; end
                for io = numel(o):-1:1
                    for ic=numel(STCONC):-1:1
                        o(io).UserData.STCONC{ic}=STCONC{ic}(o(io).idx);
                    end
                end
            end
            
            %% HEAD   requested?
            [H,varargin] = getProp(varargin,'HEAD',[]);
            if ~isempty(H)
                for io=numel(o):-1:1
                    idx_ = o(io).idx;
                    o(io).t = [H.time];
                    o(io).Dt= diff([0 o(io).t]);
                    for it=length(H):-1:1
                        o(io).UserData.H(it).values = H(it).values(idx_);
                        o(io).UserData.H(it).tstp   = H(it).tstp;
                        o(io).UserData.H(it).period = H(it).period;
                        o(io).UserData.H(it).pertim = H(it).pertim;
                        o(io).UserData.H(it).time   = H(it).time;
                        if ~isempty(o(iw).fQ)
                            o(io).UserData.H(it).value = o(io).UserData.H(it).values * o(io).fQ';
                        end
                    end
                end
            end

            %% DRAWDOWN requested?
            [DDN,varargin] = getProp(varargin,{'DD','DRAWD'},[]);
            if ~isempty(DDN)
                for io=numel(o):-1:1
                    idx_ = o(io).idx;
                    o(io).t = [DDN.time];
                    o(io).Dt= diff([0 o(io).t]);
                    for it=length(DDN):-1:1
                        o(io).UserData.DDN(it).values = DDN(it).values(idx_);
                        o(io).UserData.DDN(it).tstp   = DDN(it).tstp;
                        o(io).UserData.DDN(it).period = DDN(it).period;
                        o(io).UserData.DDN(it).pertim = DDN(it).pertim;
                        o(io).UserData.DDN(it).time   = DDN(it).time;
                        if ~isempty(o(iw).fQ)
                            o(io).UserData.DDN(it).value = o(io).UserData.DDN(it).values * o(io).fQ';
                        end

                    end
                end
            end
            
            %% Concentrations requested?
            %  Must be specified as cell array
            %  Deal with the cell array if present
            ic = strmatchi('cell',cellfun(@class,varargin,'uniformOutput',false));
            if ic(1) % cell array present
                titleStrs = varargin{ic(1)}; % varargin(ic)=[];
                i = strmatchi('head',titleStrs);
                if  i(1)
                    titleStrs(i)=[];
                    H = readDat([basename '.HDS']);
                    for io=numel(o):-1:1
                        idx_ = o(io).idx;
                        o(io).t = [H.time];
                        o(io).Dt= diff([0 o(io).t]);
                        for it=length(H):-1:1
                            o(io).UserData.H(it).values = H(it).values(idx_);
                            o(io).UserData.H(it).tstp   = H(it).tstp;
                            o(io).UserData.H(it).period = H(it).period;
                            o(io).UserData.H(it).pertim = H(it).pertim;
                            o(io).UserData.H(it).time   = H(it).time;
                            if ~isempty(o(iw).fQ)
                                o(io).UserData.H(it).value = o(io).UserData.H(it).values * o(io).fQ';
                            end

                        end
                    end
                end
                i = strmatchi({'dd','draw'},titleStrs);
                if i(1)
                    titleStrs(i)=[];
                    DDN = readData([basename,'.DDN']);
                    for io=numel(o):-1:1
                        idx_ = o(io).idx;
                        o(io).t = [DDN.time];
                        o(io).Dt= diff([0 o(io).t]);
                        for it=length(DDN):-1:1
                            o(io).UserData.DDN(it).values = DDN(it).values(idx_);
                            o(io).UserData.DDN(it).tstp   = DDN(it).tstp;
                            o(io).UserData.DDN(it).period = DDN(it).period;
                            o(io).UserData.DDN(it).pertim = DDN(it).pertim;
                            o(io).UserData.DDN(it).time   = DDN(it).time;
                            if ~isempty(o(iw).fQ)
                                o(io).UserData.DDN(it).value = o(io).UserData.DDN(it).values * o(io).fQ';
                            end
                        end
                    end
                end
                for ic=length(titleStrs):-1:1
                    C = readMT3D(sprintf('MT3D00%d.UCN',ic));
                    for io=numel(o):-1:1
                        idx_ = o(io).idx;
                        o(io).t = [C.time];
                        o(io).Dt= diff([0 o(io).t]);
                        o(io).species{ic} = titleStrs{ic};
                        o(io).NCOMP = max(o(io).NCOMP,ic);
                        for it=length(C):-1:1
                            o(io).UserData.(titleStrs{ic})(it).trpstp = C(it).trpstp;
                            o(io).UserData.(titleStrs{ic})(it).tstp   = C(it).tstp;
                            o(io).UserData.(titleStrs{ic})(it).period = C(it).period;
                            o(io).UserData.(titleStrs{ic})(it).time   = C(it).time;
                            o(io).UserData.(titleStrs{ic})(it).label  = C(it).label;
                            o(io).UserData.(titleStrs{ic})(it).tstp   = C(it).tstp;
                            o(io).UserData.(titleStrs{ic})(it).values = C(it).values(idx_);
                            if ~isempty(o(iw).fQ)
                                o(io).UserData.(titleStrs{ic})(it).value = (o(io).UserData.(titleStrs{ic})(it).values) * o(io).fQ';
                            end
                        end
                    end
                end                
            end            
            o = o([o.nr]>0);                                    
        end            

        function stamp  =stamp(o,fmt), stamp=datestr(o.created,fmt); end
        function created=get.created(o), created=o.created; end        
        function o = label(o,varargin)
            % o = observationObj.label(plotOptions); plots labels next to observations)
            for iw=1:length(o)
                text(o(iw).x,o(iw).y,sprintf('_%s',o(iw).Nr),varargin{:});
            end
        end
    end
end

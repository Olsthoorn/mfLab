
classdef wellObj
    % well= wellObj(basename, WsheetNm|well [,gr,HK [,QsheetNm [,CsheetNm],active]);
    % wellObj -- definition of class wellObj by calling
    % well = wellObj(basename,wellSheetNm) where
    %  basename is the workbook holding wellSheetNm
    %  wellSheetNm is the worksheet with a table defining the wells. See
    %  constructor by typing
    % help wellObj.wellObj.
    % wellObj defines wells for use in groundwater models. The objects hold
    % property values and share methods that facilitate using wells.
    % Some methods need a grid object. But a grid object is not necessary
    % to create well objects. This is because upon creation, well objects
    % only contain their pertinent data. Later on, when a grid object is
    % available they may be put into the grid, upon which additional grid
    % values are added to the wells, so that each well then also holds its
    % position in the grid. See the methods of wellObj and of gridObj for
    % more information.
    %
    % To immediately generate wells with concentrations and the related WEL
    % and PNTSRC array use
    % [wel,WEL,PNTSRC] = gr.well(basename,HK,sheetnmwihwells)
    %
    % SEE ALSO gridObj.well
    %
    % TO 120807
    properties
        nr;       % well nr or id (see id)
        id;       % for instance a GIS_id or equivalent of nr depending on use;
        name=''; longname='';
        x; y; z;  rw=0.25;
        ztop=0;   % ztop=top of well casing for plotting
        ix; iy; idx; iLay
        LRC       % Layer row col for the screen penetrated cells
        DZ;       % thickness of mode layers penentrated by screen
        T;        % transmissivity of layer
        fQ;       % fraction of total Q from each model layer
        Dt; t;    % stress period length and time at end of stres period
        Q; C; Cout; NCOMP; % C=specified injection conc,
        species   % names of the species/components
                  % Cout=simulated dissolved concentration in well screen Q-mixed
        parent;   % well series id to which well belongs (wellSeries(..).id)
        
        wpix=2;   % pixelwidth of wells on screen for plotting
        
        whdl=NaN(3,1);  % 3-valued handle: well casing and well screen in XSdrawing
                  %  whdl(1)=screen
                  %  whdl(2)=casing
                  %  whdl(3)=marker on casing top
        
        ITYPE=2   % well ITYPE for MT3DMS point source specfication
        
        remark = ''; % any text
        code,     % for whatever purpose
        % plotting wells, optioins:
        created   % moment at which well was create by the programm
        FaceColor = 'w';  % outline color of plotted wells in XS or YS
        EdgeColor = 'b';  % screen coor of plotted wells in XS or YS
        FaceAlpha = 1;    % transparency of wel screen
        EdgeAlpha = 1;    % transparency of lines
        marker    = 'o';        
        UserData  % for any purpose
    end
    methods
        function o=wellObj(varargin)
            % generate well objects
            % USAGE1:
            %       wel = wellObj(nr,x,y,z[,rw]); % basic way of generating wells (one by one)
            % USAGE2:
            %       well= wellObj(basename, well [,gr,HK [,QsheetNm [,CsheetNm[,use]]);
            %
            %       This call serves to update existing wells (well) with new grid, Q and conc. info
            %
            % USAGE3:
            %       well= wellObj(basename, wellSheetNm [,gr,HK [,Qloc [,Cloc]]);
            %
            %       This call generates wells using basename,wellSheetNm, and
            %             with optional arguments 3 (grid) and 4 (HK | TRAN)
            %             puts them into the grid.
            %             With further arguments it fetches flow
            %             and concentration data from a worksheet.
            %             Qloc is a string or a cell array with strings
            %             that tells the location of the well flows.
            %             If complete it consists of the name of the
            %             worksheet containing the data in a table with
            %             values for all stress periods vertically. Plus
            %             the prefix of the columns with the flow data of
            %             the wells. The default worksheet is 'PER' and the
            %             default header is Q or Q_. You must specify
            %             something, so generally use 'PER' for it. The
            %             data for the wells are then retrieved from this
            %             spreadsheet under column headers Q1, Q2, Q3 etc.
            %             where the number refers to the well nr, i.e.
            %                  (well(iw).nr)
            %             examples {'PER','Qaxi'}
            %             If there is only 1 header that matches the prefix
            %             exactly, all wells will get the data from that
            %             column.
            %             The Cloc arguments is similar. It serves to fetch
            %             concentration data for the wells. The default
            %             spreadsheet is, again, 'PER' and the default
            %             header prefix is C or C_. Ideally use the species
            %             name as the header. Then you can get the data
            %             like this:
            %             {'PER','Salinity','Temp'}
            %             In that case the salinity of well 5 will be
            %             fetched from column 'Salinity5' in the 'PER'
            %             worksheet, and the temmperature from column
            %             'Temp5' from the same worksheet.
            %
            %             To only select a subset of the wells defined in
            %             the spreadsheet wellSheetNm or given as wellObj
            %             under USAGE2, you can use the string,value pair
            %             'active',values  or 'use',values in the command
            %             line at any location. values are either logicals
            %             identifying which wells are active, or numbers
            %             that refer to the well numbers or a character
            %             expression that is applied on the name of each
            %             well through regexp() funcion (regular
            %             expression).
            %             instead of the string value pair you can specify
            %             activeWells as the 6th argument.
            %
            %       some shortcuts, however, wellObj is preferred because of
            %           better error checking of the input
            %
            %       well = well.inGrid(gr,HK);                % puts wells into the grid
            %       well = well.getQ(basename[,WsheetNm]);    % extracts Q from basename sheet PER
            %       well = well.getCin(basename[,WsheeteNm]); % extracts Cin from basename, sheet PER
            %       well = well.getQCin(basename[,WsheetNm]); % extracts both Q and Cin from sheet PER
            %
            %       EXAMPLES:
            %       wellObj(basename,'wells');
            %       wellObj(basename,'wells',gr,HK,'PER','PER','use',[1 2 3]);
            %       wellObj(basename,'wells',gr,HK,'Qaxi','PER','use',[1 2 5 7]);
            %       wellObj(basename,'well,'active','Am*',{'PER','Qflat'},{'Temp','Salinity'})
            % SEE ALSO: gridObj.well (shortcut)
            %
            % TO 120807
            
            if nargin<2, return; end

            % fill the well table if cells are empty or use NaNs ?
            [fill,varargin]        = getProp(varargin,'fill'  ,true);
            [indexFld,varargin]    = getProp(varargin,'index','');
            index = ~isempty(indexFld);
            
            [activeWells,varargin] = getProp(varargin,{'use','active'},[]);
            if numel(varargin)>6, activeWells = varargin{7}; varargin(7:end)=[]; end
            
            [gr,varargin] = getType(varargin,'gridObj',[]); 

            % USAGE1:
            if nargin == 6 && all(cellfun(@isnumeric,varargin))

                for j=numel(varargin{1}):-1:1
                    o(j).nr   = varargin{1};
                    o(j).id   = o(j).nr;
                    o(j).x    = varargin{2}(j);
                    o(j).y    = varargin{3}(j);
                    o(j).z(1) = varargin{4}(j);
                    o(j).z(2) = varargin{5}(j);
                    o(j).rw   = varargin{6}(j);
                    o(j).name =sprintf('Well%03d',o(j).nr);
                    o(j).created=now;
                end
                return
            end

            
            [basename,varargin] = getNext(varargin,'char',[]);
            if isempty(basename)
                error('%s: First argument must be the basename (type char, not <<%s>>.',...
                    mfilename,class(varargin{1}));
            end

            [wellSheetNm,varargin] = getNext(varargin,'char',[]);

            if ~isempty(wellSheetNm)
                % USAGE2;  % well sheet name given
                
                %% get the data from the spreadsheet with the given name
                members = fieldnames(o);

                use = [];
                table = getTable(basename,wellSheetNm,members,'Horizontal',use,'fill',fill);
                field = fieldnames(table);

                for j=1:numel(field)
                    for iw=numel(table):-1:1
                        o(iw,1).(field{j}) = table(iw).(field{j});
                    end
                end

                % deal with well radius separately
                if isfield(table,'UserData')
                    UDfield = fieldnames(table(1).UserData);  % UDfield = User Data field names
                    ic = strmatchi('diam',UDfield);
                    if ic(1)
                        for iw=length(o):-1:1
                            o(iw).rw = o(iw).UserData.(UDfield{ic(1)})/2;
                        end
                    end
                    ic = strmatchi('Nr',UDfield,'exact');
                    if ic(1)
                        for iw=length(o):-1:1
                            o(iw).nr = o(iw).UserData.(UDfield{ic(1)});
                        end
                    end
                end
                
                %% deal with possible fields like z1 and z2 separately,
                % zScreen..., that all join into screen z coordinates.
                % TODO
                %   This is meant for curved, i.e. multi section screens and may need some work and
                %   precision to be robust
                if isempty(o(iw).z)
                    L = regexpis(UDfield,'z[1-9]|zscr'); % indices of fields z1 z2 z3 z4 ... etc
                    if isempty(L)
                        error('%s: Add columns z1 z2 to well specification table to specify top and bottom of well screen',mfilename);
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

            else % USAGE3; % well objects given

                [o, varargin] = getNext(varargin,{'wellObj','MNW1Obj','MNW2Obj'},[]);
                if isempty(o)
                    error(['%s: No wells found, make sure the second arg is either\n',...
                        '1) the name of the spreadsheet with the well info\n',...
                        '2) a well object array'],mfilename);
                end
            end
            
            %% Remove wells with zero or negative well numbers (or NaNs)
            o = o([o.nr]>0);
            
                        %% ======= Only keep active wells =============================
            if ~isempty(activeWells)
                if islogical(activeWells)
                    o = o(activeWells);
                elseif isnumeric(activeWells)
                    o = o(ismember([o.nr],activeWells));
                elseif ischar(activeWells)
                    for io=numel(o):-1:1
                        if ~regexpi(o(io).name,activeWells,'once')
                            o(io)=[];
                        end
                    end
                else
                    % ignore use
                end
            end
            
            %% Any wells left ?
            if isempty(o)
                error(['%s: no wells left: possible reasons:\n',...
                    '1) There is no column with header ''nr'' in the well specification table\n',...
                    '2) All fields in column ''nr'' are empty\n',...
                    '3) All fields in column ''nr'' are <=0\n',...
                    ' because well with nr<=0 or empty are eliminated!\n',...
                    'REMEDY: Verify the above points.'],...
                    mfilename);
            end
            
            %% Fill in well.id if still empty
            for iw = numel(o):-1:1
                if isempty(o(iw).id)      , o(iw).id       = o(iw).nr;   end
                if isempty(o(iw).name)    , o(iw).name     = sprintf('well%d',o(iw).nr); end
                if isempty(o(iw).longname), o(iw).longname = o(iw).name; end
            end

            % At this pont get HK from input line if it exists
            [HK,varargin] = getNext(varargin,'double',[]);

            %% ======= Get flows from QsheetNm ============================
            [QsheetNm_or_QCol,varargin] = getNext(varargin,{'char','cell'},[]);
            if ~isempty(QsheetNm_or_QCol)
                if ~index
                    o = o.getQ(basename,QsheetNm_or_QCol,'fill',fill);
                else
                    o = o.getQ(basename,QsheetNm_or_QCol,'index',indexFld,'fill',fill);
                end
            end
                        
            %% ======= Get input concentrations from CsheetNm =============
            [CsheetNm_or_Ccol,varargin] = getNext(varargin,{'char','cell'},[]);
            if ~isempty(CsheetNm_or_Ccol)
                if ~index
                    o = o.getCin(basename,CsheetNm_or_Ccol,'fill',fill);
                else
                    o = o.getCin(basename,CsheetNm_or_Ccol,'index',indexFld,'fill',fill);
                end
            end

            if ~isempty(varargin)
                fprintf('%s: varargin variables not all used',mfilename);
            end
            
            %% =======  Continue: Put wells into grid =====================
            if isempty(HK) || isempty(gr), return; end

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
            
            if ~isa(o,'MNW1Obj') && ~isa(o,'MNW2Obj') && numel([o.nr])~=numel(unique([o.nr]))
                msgId      = 'mfLab:wells:notUnique';
                warning('on',msgId);
                beep();
                warning(msgId,['One or more wells are not unqiue, as they have the same well number!\n',...
                    'Notice that their screens will be combined. This may not be what you indented.']);
                warning('off',msgId);
            end

            o = o.toGrid(gr,HK);  % toGrid(gr,HK) (original HK)
            
        end
            
        function X   = X(o)
            % get first x coordinate of wells
            for io=numel(o):-1:1
                X(io) = o(io).x(1);
            end
        end
        function Y = Y(o)
            % get first y coordinate of wells
            for io=numel(o):-1:1
                Y(io) = o(io).y(1);
            end
        end
        function Z_top = Z_top(o)
            % get highest z-coordinate of well
            for io=numel(o):-1:1
                Z_top(io)= o(io).z(1);
            end
        end
        function Z_bot = Z_bot(o)
            % get lowest z-coordinate of well
            for io=numel(o):-1:1
                Z_bot(io)=o(io).z(end);
            end
        end
        function Idx = Idx(o)
            % Get first index of every wellObj in a vector
            % USAGE: Idx= wellObj.Idx;
            for io=numel(o):-1:1
                Idx(io)=o(io).idx(1);
            end
        end
        
        function Ix= Ix(o)
            for io=numel(o):-1:1
                Ix(io) = o(io).ix(1);
            end
        end
        
        function Iy= Iy(o)
            for io=numel(o):-1:1
                Iy(io) = o(io).iy(1);
            end
        end
        
        function stamp  =stamp(o,fmt), stamp=datestr(o.created,fmt); end
        function created=get.created(o), created=o.created; end
                
        function Cout=get.Cout(o), if isempty(o.Cout), Cout=NaN(size(o.C)); else Cout=o.Cout; end; end
        function C   =get.C(o),    if isempty(o.C),    C   =NaN(size(o.Q)); else C   =o.C;    end; end
        
%         function MRout = get.MRout(o), MRout= o.Q.*o.Cout; end    % mass rate out of well
%         function MRin  = get.MRin (o), MRin = o.Q.*o.C   ; end    % mass rate in  of well
%         function Mout  = get.Mout (o), Mout = o.MRout.*o.Dt; end % mass rate multiplied by length of SP
%         function Min   = get.Min  (o), Min  = o.MRin .*o.Dt; end % same
        
        function leg = plotQ(   o      ,varargin)
            % wellObj.plotQ(ax,plotOptions); plots well flow Q accoring to
            % plotOptions: eg.
            % well(10).plotQ('lineStyle','--','color','r');

            [ax         ,varargin] = getNext(varargin,'axis',gca);

            set(ax,'nextplot','add');

            set(get(ax,'title'),'string','Q')

            leg{1,length(o)} = '';
            for iw=1:length(o)
                if isempty(varargin)
                    plot(o(iw).t,o(iw).Q  ,mf_color(iw));
                else
                    plot(o(iw).t,o(iw).Q  ,varargin{:});
                end
                leg{iw} = o(iw).name;
            end
        end

        function leg = plotC(o,varargin)
            % wellObj.PlotC(ax,speciesName,iComp,varargin); plot well input concentration of
            % component iComp according to plotOptions.

            [ax         ,varargin] = getNext(varargin,'axis',gca);
            [speciesName,varargin] = getNext(varargin,'char','');
            [iComp      ,varargin] = getNext(varargin,'double',1);

            if isempty(speciesName)
                speciesName = sprintf('speciesName%d',iComp);
            end

            set(ax,'nextplot','add');

            set(get(ax,'title'),'string',['Cfix ',speciesName])

            leg{1,length(o)} = '';
            for iw=1:length(o)
                if isempty(varargin)
                    plot(ax,o(iw).t,o(iw).C(iComp,:),mf_color(iw));
                else
                    plot(ax,o(iw).t,o(iw).C(iComp,:),varargin{:});
                end
                leg{iw} = o(iw).name;
            end
        end
        
        function leg = plotCout(o,varargin)
            %PLOTCOUT -- plots output concentration of well
            % for component iComp and according to plotOptions.
            %
            % USAGE: 
            %    leg = wellObj.plotCout(ax[,speciesName],iComp,plotOptions);
            %
            % TO 120101
            
            [ax         ,varargin] = getNext(varargin,'axis',gca);
            [speciesName,varargin] = getNext(varargin,'char','');
            [iComp      ,varargin] = getNext(varargin,'double',1);

            if isempty(speciesName)
                speciesName = sprintf('speciesName%d',iComp);
            end

            set(ax,'nextplot','add');

            if isempty(get(get(ax,'title'),'string'))
                set(get(ax,'title'),'string',['Cout ',speciesName]);
            end

            leg=[];
            for iw=1:length(o)
                if isempty(varargin)
                    plot(ax,o(iw).t,o(iw).Cout(iComp,:),mf_color(iw));
                else
                    plot(ax,o(iw).t,o(iw).Cout(iComp,:),varargin{:});
                end
                leg{iw} = o(iw).name;
            end
            legend(leg{:},2);
        end
        
        function plotCin( o,varargin)
            % wellObj.plotCin(ax,iComp,plotOptions); plots input
            % concentration of component iComp and according to plotOptions
            % 'color','k' etc.
            
            [ax,    varargin] = getNext(varargin,'axis',gca);            
            [iComp, varargin] = getNext(varargin,'double',1);

            set(ax,'nextplot','add');
            
            for iw=1:length(o)
                if isempty(varargin)
                    plot(o(iw).t,o(iw).C(iComp,:),mf_color(iw));
                else
                    plot(o(iw).t,o(iw).C(iComp,:),varargin{:});
                end
            end
        end
        
        function o = plotXY(o,varargin)
            % well = well.plotXY([ax|it],plotOptions); plots well locations according to
            % plotoptions

            [ax,varargin] = getNext(varargin,'axis',[]);
            [it,varargin] = getNext(varargin,'double',1);

            if isempty(ax)
                if it==1
                    ax = gca;
                end
            else
                it = 1;
            end
                                 
            if isempty(varargin), varargin={'visible','on'}; end

            if it==1;
                set(ax,'nextplot','add');
                lineSpecGiven = false;
                if numel(varargin)>0 && isLineSpec(varargin{1})
                         lineSpec = varargin{1};
                         varargin(1)=[];
                         lineSpecGiven = true;
                end
                for iw=1:length(o)
                    o(iw).whdl(3) = NaN;
                    if lineSpecGiven
                        o(iw).whdl(3)=plot(ax,o(iw).x,o(iw).y,lineSpec,varargin{:});
                    else
                        o(iw).whdl(3)=plot(ax,o(iw).x,o(iw).y,...
                            'markerFaceColor',o(iw).FaceColor,...
                            'markerEdgeColor',o(iw).EdgeColor,varargin{:});
                    end
                end
            end
                        
            for iw=1:numel(o)
                if numel(o(iw).Q)>=it
                    if isnan(o(iw).Q(it))
                        set(o(iw).whdl(3),'visible','off');
                    else
                        set(o(iw).whdl(3),'visible','on');
                    end
                    if o(iw).Q(it)==0
                        set(o(iw).whdl(3),'markerFaceColor',grey);
                    else
                        set(o(iw).whdl(3),'markerFaceColor',o(iw).FaceColor);
                    end
                else
                    %set(o(iw).whdl(3),'visible','off');
                end
            end
        end
        function o = label(o,varargin)
            % well = wellObj.label(plotOptions); plots labels next to wells)
            for iw=1:length(o)
                text(o(iw).x,o(iw).y,sprintf('_%s',o(iw).Nr),varargin{:});
            end
        end
        function o = plot3D(o,varargin)
            % well = wellObj.plot3D([ax,]plotOptions); plots wells in 3D according to
            % plotOptions.
            
            [ax,varargin] = getNext(varargin,'axis',[]);
            if isempty(ax)
                ax = gca;
            else
                for iw=1:numel(o)
                    o(iw).whdl =NaN(3,1);
                    if isempty(o(iw).EdgeColor), o(iw).EdgeColor = grey; end
                    if isempty(o(iw).FaceColor), o(iw).FaceColor = grey; end
                end
            end
            
            if isempty(varargin)
                varargin = {'visible' 'on'};
            else
                [LS,c,m,L] = isLineSpec(varargin{1});
                if LS
                    if c
                        varargin = [varargin, 'color', c];
                    end
                    if m
                        varargin = [varargin, 'marker', m];
                    end
                    if ~isempty(L)
                        varargin = [varargin, 'lineStyle', L{:}];
                    end                   
                    varargin = varargin(2:end);
                end
            end

            o = o.plotScreen3D( ax,varargin{:},'lineWidth',3);
            o = o.plotCasing3D( ax,varargin{:},'lineWidth',0.1);
            o = o.plotScreenTop(ax,varargin{:});            
        end
        function o = plotScreen3D(o,varargin)
            % well = wellObj.Screen3D(ax|it,plotOptions); plots wells screens in 3D according to
            % plotOptions.
            
            [ax,varargin] = getType(varargin,'axis',gca);
               
            if isempty(varargin)
                varargin = {'visible','on'};
            end
                
            if isnumeric(varargin{1})
                it = varargin{1};
                varargin(1)=[];
            else
                it=1;
            end
            
            if ~strmatchi('lineWidth',varargin)
                varargin = [varargin, {'lineWidth' 2}];
            end

            if it==1
                for iw=1:numel(o)
                    o(iw).whdl(1) = plot3(ax,[o(iw).x(1) o(iw).x(1)],[o(iw).y(1) o(iw).y(1)],o(iw).z([1 end]),varargin{:});
                    o(iw).EdgeColor = get(o(iw).whdl(1),'color');
                end
            end
            
            for iw=1:numel(o)
                % whdl(1) = screen
                if isnan(o(iw).Q(it))
                    set(o(iw).whdl(1),'visible','off');
                else
                    set(o(iw).whdl(1),'visible','on');
                end
                if o(iw).Q(it) ==0
                    set(o(iw).whdl(1),'color',grey,varargin{:});
                end
            end
        end
        function o = plotCasing3D(o,varargin)
            % well = wellObj.Casing3D(ax|it,plotOptions); plots wells casings in 3D according to
            % plotOptions.
            
            [ax,varargin] = getType(varargin,'axis',gca);

            if isempty(varargin)
                varargin = {'visible' 'on'};
            end
            
            if isnumeric(varargin{1})
                it = varargin{1};
                varargin(1)=[];
            else
                it=1;
            end
                        
            if it==1
                for iw=1:numel(o)
                    o(iw).whdl(2)= plot3(ax,[o(iw).x(1) o(iw).x(1)],[o(iw).y(1) o(iw).y(1)],[o(iw).ztop o(iw).z(1)],varargin{:});
                    o(iw).EdgeColor = get(o(iw).whdl(2),'color');
                end
            end
            
            for iw=1:numel(o)
                % whdl(2) is casing, whdl(3) is marker on top of casing
                if isnan(o(iw).Q(it))
                    set(o(iw).whdl(2),'visible','off');
                else
                    set(o(iw).whdl(2),'visible','on');
                end
            end
        end
        function o = plotScreenTop(o,varargin)
            % well = wellObj.plotScreenTop(ax|it,plotOptions); plots wells tops in 3D according to
            % plotOptions.
            
            [ax,varargin] = getType(varargin,'axis',gca);

            if isempty(varargin)
                varargin = {'visible' 'on'};
            end
            
            if isnumeric(varargin{1})
                it = varargin{1};
                varargin(1)=[];
            else
                it=1;
            end
            
            [color          ,varargin] = getProp(varargin,'color',[]);                
            [markerEdgeColor,varargin] = getProp(varargin,'edgeColor',color);
            [markerFaceColor,varargin] = getProp(varargin,'faceColor',color);
            [Marker         ,varargin] = getProp(varargin,'marker'   ,'o');
            
            for iw = numel(o):-1:1
                if isempty(markerFaceColor)
                    if ~isempty(color)
                        o(iw).FaceColor = color;
                    end
                else
                    o(iw).FaceColor = markerFaceColor;
                end
                if isempty(markerEdgeColor)
                    if ~isempty(color)
                        o(iw).EdgeColor = color;
                    end
                else
                    o(iw).EdgeColor = markerEdgeColor;
                end
                if ~isempty(Marker)
                    o(iw).marker = Marker;
                end
            end
            
            if isaxis(ax)
                for iw=1:numel(o)                                        
                    o(iw).whdl(3)= plot3(ax, o(iw).x(1),...
                                             o(iw).y(1),...
                                             o(iw).ztop+100*eps,...
                                             varargin{:},...
                                             'marker'         ,o(iw).marker,...
                                             'markerEdgeColor',o(iw).EdgeColor,...
                                             'markerFaceColor',o(iw).FaceColor);
                end
            end
            
            for iw=1:numel(o)
                % whdl(2) is casing, whdl(3) is marker on top of casing
                if isnan(o(iw).Q(it))
                    set(o(iw).whdl(3),'visible','off');
                else
                    set(o(iw).whdl(3),'visible','on');
                end
                if o(iw).Q(it) ==0
                    set(o(iw).whdl(3),varargin{:},'markerFaceColor',grey);
                else
                    set(o(iw).whdl(3),varargin{:},'markerFaceColor',o(iw).FaceColor);
                end
            end
        end

        function write(o,fid,iPer)
           %% well.write(fid,iPer) -- write this well to file for this stress period
           %  using in writeBCN.m called by mf_setup
           %  TO 120629
           for iw=1:length(o)
               if ~isnan(o(iw).Q(iPer)) && o(iw).Q(iPer)~=0,
                   for iCell = 1:length(o(iw).idx)
                       fprintf(fid,'%10d%10d%10d%10g\n',...
                           o(iw).LRC(iCell,:),...
                           o(iw).fQ(iCell)*o(iw).Q(iPer).');
                   end
               end
           end
        end
       function writePNTSRC(o,fid,iPer)
           %% WEL.writePNTSRC(fid,iPer) -- write PNTSRC for this WEL to file for this stress period
           %  using in writeBCN.m called by mf_setup
           %  TO 120629
           fmt = sprintf('%%10d%%10d%%10d%%14g%%3d%s\n',repmat('%14g',[1,size(o.C,1)]));
           for iCell=1:length(o.idx)
               fprintf(fid,fmt,...
                   o.LRC(iCell,:),...
                   o.C(1,iPer),...
                   o.ITYPE,...
                   o.C(:,iPer).');
           end
       end
       
       function zoneArray = ZONE(o,zoneArray)
           % zoneArray = well.ZONE(zoneArray) -- add well(iw).nr to zoneArray
           for iw = 1:numel(o)
               zoneArray(o(iw).idx) = o(iw).nr;
           end
       end
       
       function o = removeUnderscores(o)
           % remove underscores from names of wells
           % capitalize the letter after the underscore.
           for iw=numel(o):-1:1
               if ~isempty(o(iw).name)
                   I = o(iw).name == '_';                    
                   if any(I),
                       J = I(I<numel(o(iw).name))+1;
                       o(iw).name(J)=upper(o(iw).name(J));
                       o(iw).name(I)=[];
                   end
               end
               if ~isempty(o(iw).longname)
                   I = o(iw).longname=='_';
                   if any(I)
                       J = I(I<numel(o(iw).name))+1;
                       o(iw).longname(J)=upper(o(iw).longname(J));
                       o(iw).longname(I)=[];
                   end
               end
           end
       end

    end
end

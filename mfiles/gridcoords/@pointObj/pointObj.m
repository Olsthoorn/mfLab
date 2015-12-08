classdef pointObj
            %POINTOBJ: or lineObj or area2Obj generates point objects, line objects
            % or area2 objects defined in a spreadsheet to define boundary conditions
            % or stresses like rivers, drains, general head boundaries and
            % boundaries with specified head and flux.
            %
            % pointObj is a superclass of lineObj which is a superclass of area2Obj
            %          they are consistent.
            %
            % SEE ALSO wellObj pointObj lineObj area2Obj
            %
            % Usage instructions for pointObj are also valid for lineObj
            % and area2Obj. The help instructions for the lineObj and
            % area2Obj refer to those of the pointObj given here.
            %
            % USAGE:
            %       point= pointObj([gr,] basename,sheetNm, type[,colName toSelect],...
            %         ['active',perSheetColName,] ['conc',{speciesName1 speciesName2 ...}],...
            %         [value1|plane1,[value1|plane2[,...]]])
            %
            % Alternative: with basename='' and sheetNm='' and a struct
            % given as a parameter pair value and rest unchanged to allow
            % specifying coordinate data and attributes through a struct
            % instead of through Excel. This allows for instance use of
            % shapes directly and can be handy in other circumstances as
            % well.
            %       point= pointObj([gr,] '', '',type, ... 'struct',struct, ...);
            %
            %     See for use of 'struct',struct instead of
            %     basename,sheetNm and Excel the instructions as the very
            %     end of this help info.
            %
            % Notice: instead of pointObj one may also use lineObj and
            %         area2Obj. All these objects can be used to specify
            %         stress boundaries for MODFLOW and to specify
            %         stress concentrations for MT3DMS and SEAWAT.
            %         pointObj, lineObj, area2Obj are completely defined by
            %         their vertices defined in a spreadsheet table as
            %         explained below. This makes for efficient and
            %         consistent input for the different object types.
            %         A pointObj consists of a single vertix (one
            %         spreadsheet line)
            %         A lineObj and area2Obj are difined by at least two or
            %         three vertices (lines) in the data spreadsheet.
            %         With a lineObj the vertices define a river etc.
            %         With an area2Obj the vertices enclose an area. This
            %         area will automatically be closed if the beginning
            %         and end vertices differ.
            %         So you can freely exchange pointObj for lineObj and
            %         areaObj below.
            %
            % The most basic call would be
            %       point = pointObj(basename,sheetNm);
            %       point = pointObj(basename,sheetNm,type)
            %
            % Without a gridObj specified, the data are read from Excel
            % workbook basename, sheet sheetNm, but no information about the
            % position of the object within the model grid will and can be
            % generated. So methods that require that information will
            % fail. But the object can plot itself, e.g.:
            %
            %       point.plot('ro','markerSize',8)
            %       point.plot(ax,'ro','markerSize',8)
            %
            % INPUT:
            %   basename    is the name of the Excel workbook to read the data from.
            %   sheetNm     is the name of the sheet in the workbook holding the data.
            %     The data in the spreadsheet consist of a table where each line
            %     of the table specifies a pointObj or the vertex of a lineObj
            %     or a vertex of an area2Obj.
            %     A lineObj or an area2Obj, therefore, is defined by
            %     their vertices, e.g. a set of lines of the table in the sheet.
            %   type       is stress type of the object(s). This type must be
            %     recognized as one of
            %          {'NIL' 'WEL'|'FLUX' 'DRN' 'GHB' 'RIV' 'CHD'}
            %    'WEL' and 'FLUX' are equivalent. The set of implemented types
            %     may be extended in the future.
            %     NIL is a default if no stress type is given.
            %
            % The rest are optional arguments:
            % colName,toSelect    is a pair of strings denoting a column name
            %    in the data and the  values to select. toSelect should be the
            %    first letters of a string common to the values used to
            %    extract the vertices from the spreadsheet. If this argument
            %    pair is omitted, all lines will be extracted. They will be
            %    split into distinct objects based on the column 'name' which
            %    is obligatory. The values in this column must be either
            %    numeric or character strings
            %  Example:                                %
            %       point = pointObj(gr,basename,sheetNm,type,'region','dunes');
            %
            % Only the lines will be read that have in teh column with header
            % 'region'  lines with stings that commence with 'dune' hence
            % 'dune*'.
            % These lines will subsequently be collated into disctinct
            % objects based on the uniqueness in the string values in the
            % specified column.
            % This would include strings like 'dunes1','dunes2','dunesZandvoort'
            % etc. The unique strings in resulting column 'region' define the
            % distinct resulting objects. For instance if we have 8 lines
            % with string 'dunes2', these lines will be interpreted as subsequent
            % vertices of the object 'dunes2'.
            % Of course, for a pointObj each line should be unique as there can be
            % only one vertex per objects, but for lineObj there has to be at least
            % 2 vertices per object and for area2Obj at least 3.
            %
            % To indicate when an objects is active or incative, add the
            % following arguments to the intput list
            %   ...,'active',perSheetHdrName,...
            % Upon generating input for MODFLOW etc, the column in the PER
            % sheet with header equal to the specified perSheetHdrName will
            % then be used to determine in which of the stress periods the
            % object(s) is(are) active. A 0 will be inactive, and any other
            % value active.
            % To inciated the opposite preprend perSheetHdrName with a '~':
            %  ...,'active','~working',...    versus
            %  ...,'active','working',...
            %
            %  oldSea   = pointObj(gr,basename,sheetNm,'CHD','region','newDunes',...
            %                'active','~realized',...)
            %  newDunes = pointObj(gr,basename,sheetNm,'DRN','region','newDunes',...
            %                'active','realized',...)
            %
            % When running MT3DMS or SEAWAT the concentrations of the
            % stress points also need to be given to specify what
            % concentration the program should use when water enters the
            % model at the stress points. This is done by adding the string
            % 'conc' to the intput line followed by the names of the
            % species like  ...,'conc',{'salinity','temp',...},....
            % braces { } are only required when more than one species are
            % simulated. These names refer to columns in the spreadsheet
            % providing the values for each of the vertices in case that the
            % columes are numeric. When the columns indicated by the
            % concentration headers contain strings, only the firs string
            % is used. It is interpreted as the header of a column in the
            % PER worksheet to get values from that are distinct per stress
            % period.
            %  oldSea   = pointObj(gr,basename,sheetNm,'CHD','region','newDunes',...
            %         'active','~realized','conc',{'TDS','temp'},...)
            %  newDunes = pointObj(gr,basename,sheetNm,'DRN','region','newDunes',...
            %         'active','realized','conc',{'TDS','temp'}...)
            %
            % one may also directly specify numeric values for concentrations as
            %
            %  newDunes = pointObj(gr,basename,sheetNm,'DRN','region','newDunes',...
            %         'active','realized','conc',{'TDS=0.001','temp=15'}...)            
            %
            % EXAMPLES:
            %       piez= pointObj(basename,sheetName,gr,'CHD',DEM,DEM)
            %       riv = pointObj(basename,sheetName,gr,'RIV',DEM,[],DEM)
            %       drn = pointObj(basename,sheetName,gr,'DRN',DEM}
            %       ghb = pointObj(basename,sheetName,gr,'GHB',DEM}
            %       chd = pointObj(basename,sheetName,gr,'CHD',DEM,DEM)
            %       well= pointObj(basename,sheetName,gr,'WEL',MULTARRAY)
            %       pnts= pointObj(basename,sheetName,gr);
            %
            % The table in the spreadssheet has a header with strings
            % indicating the contents in the column below each header.
            % The number of columns is arbitrary, but some columns are obligatory
            %   -- columns defining the vertex coordinates i.e.
            %         'Lat' 'Lon' and or 'x' 'y'
            %  Lat Lon will be automatically converted to UTM x and y if
            %  both columns 'x' and 'y' are missing.
            %  --- either a numeric column 'Id' or a string column 'name'
            %  are required. Instead column 'name' can be replaced by
            %  another column as explained below.
            %  --- depending on the stress type, one to 3 extra columns are
            %  obligatory, as they will be used to write out dat for
            %  MODFLOW etc. For instance, type 'WEL' requires one column
            %  with the flow values. So a column with header 'WEL','FLUX'
            %  or 'Q' is required. For instance, type 'CHD' requires the
            %  head at the beginning and the end of each stress period.
            %  Therefore, two columns are required with headers
            %  respectively 'H1|Hstart|H' and 'H2|Hend|H'. If only a column
            %  with header H is found, then it will be used for both H1 and
            %  H2. Notice that if the respective column is numeric, that
            %  value will be used for all stress periods. But if it is a
            %  string, this string will be interpreted as the header in the
            %  PER worksheet and the corresponding stress period values
            %  will be used for the distinct stress periods. This is the
            %  general interpertation. Stress type 'RIV' requires three
            %  inputs, the river stage, the conductance and the river
            %  bottom elevation. Therefore three columns must be present to
            %  supply these values, resp: 'H|stage|', 'C', 'Hb*' where '*'
            %  can be zero or more characters. The the interpreation of
            %  use of numeric or string values are as explained before.
            %  Likewise, 'type' 'DRN' requires 'H|elevation' and 'C'
            %            'type' 'GHB' requires 'H|stage'     and 'C'
            %  Notice that the respective columns should be either numeric
            %  or of character type, not mixed.
            %
            %
            % where the optional arguments are a legal matlab linSpec
            % followed by any plot options that are legal to the basic
            % matlab plot function.
            % ax may be supplied, it must be a legal axis handle. Default
            % is gca.
            % For other plot options see the methods of the pointObj
            % lineObj and area2Obj.
            %
            % To generate grid information, add a a variable of class gridObj
            % to the input arguments. It is here called gr.
            % The position of this gridObj in the list of arguments
            % is immaterial, as it will be recognized automatically amongst
            % all other arugments:
            %       point = pointObj(gr,basename,sheetNm);
            %       point = pointObj(basename,sheetNm,type,gr);
            %
            %
            % The position of gr in the input is immaterial, however, if
            % the gridObj is omitted, the resulting lineOb is not put in the
            % grid, i.e. has no information about the model grid.
            % To put lines into a (new) grid
            %       riv = riv(gr);
            %       drn = drn(gr);
            %       ghb = ghb(gr);
            %       chd = chd(gr);
            %       flux=flux(gr);
            %
            %
            % 'id|name','x|Lon','y|Rel','z|zRel' that define the vertices of
            % the polyline defining the point or area object.
            % And some columns are required by the type of stress:
            % You can choose an arbitray column to be interpreted as y
            % coordinates by using the name property pair 'y',selectHdr inthe
            % input, where selectHdr is the name of the column in the input
            % spreadsheet holding the y-values. This is useful when working
            % in cross sections. One could use 'y','z' pair in the input to
            % interpret the column with header 'z' as y coordinates in the
            % object.
            %
            % Grid verices have an elevation. This elevation can be
            % specified either in abolute elevation values or in relative
            % elevation values. To define vertices elevations absolutely,
            % make sure there is a column with header 'z' in the
            % spreadsheet. To define vertices with a relative elevation
            % either specify a column with header 'zRel' or no column 'z' and
            % 'zRel' at all.
            % The values in column 'zRel' will be interpreted as relative
            % to the top of the model grid. Leaving out both the column 'z'
            % and 'zRel', will be intepreted as zRel=0, i.e. vertices at
            % ground surface.
            % Notice that when using absolute values combined with vertices
            % that are apart farther than one cell, there is a chance|risk
            % that part of the line spanning the vertices is outside the
            % model (e.g. above ground surface). This will result in only
            % part of the boundary line taken into account, with that part
            % outside the model ignored. To prvent this, use relative
            % z-coordinates. This will ensure that the spatial
            % interpolation is done within a grid with regular (relative)
            % z-coordinates, i.e. where all z-columns are equal.
            %
            % Dimensions of C (concuctance).
            % The dimension of the specified value of the vertex
            % conductance depends on the object class. For pointObj the
            % dimension is m2/d, for lineObj it is m/s and for area2Obj it
            % is 1/d. That is, if the conductance is multiplied by a head
            % difference, one obtains for pointObj the total flow for the cell in
            % which the vertex resides, for lineObj one obtains the flow
            % per m length of the lineObj intersecting the cells (m2/d), and for
            % area2Obj the flow per unit area (m/d). It is the
            % responsibilit of the user to make sure that the dimensions
            % are consistent.
            %
            % Notice that values in cells that are intersected by line
            % objects are interpolated from the supplied vertex values.
            %
            % For area2Obj, first the cells on the contour of the area are
            % interpolated from its vertices, after making sure that the
            % area is closed. Then, in a second step, the values interior
            % to the area2Obj are interpolated from its contour values by a
            % Laplacian interpolator. (i.e. as if it is groundwater with
            % its contour as fixed head).
            %
            % Extra arguments may be added to the input that allow a spatial
            % interpretation of the columns used for generating the stress
            % inputs for MODFLOW, MT3DMS and SEAWAT.
            % These inputs are either empty, scalars or are numeric arrays
            % having the size of a model layer. The number
            % of extra agruments is equal to the number of arguments
            % required by the different stresses: i.e. 1 for 'WEL|FLUX', 2 for
            % 'DRN|GHB', 3 for 'RIV' and 2 for 'CHC'.
            % These arguments are, therefore planes that affect the vertex
            % values in the columns of the spreadsheet as follows:
            % If the column to which a plane is attributed is a head or an
            % elevation, the head or elevation in the spreadsheet columns
            % are added to the plane. This makes them relative to the plane.
            % For instance, if DEM is given as a plane, the values for the
            % elevation 'H' in the spreadsheet will be relative to this plane
            % and computed. A plane that is attributed to columns that are
            % neither an elevation nor a head, i.e. only the
            % conductance column (C) or the column WEL|FLUX|q, than the plane
            % values are considered a spatial multiplier. This allows spatially varying
            % time-dependent data. E.g.
            %
            % area1 = areaObj(gr,basename,sheetNm,'FLUX', multiPlierArray);
            % area2 = areaObj(gr,basename,sheetNm,'RIV',  DEM,[],DEM)
            %
            % In the first example, the values in the column 'Q|FLUX|WEL'
            % in the spreadsheet will be multiplie by the values in the
            % multiplierArray (size gr.Ny,gr.Nx) and the values that fall
            % inside the areaObj will be attributed to it. This allows
            % distributing a given flux value in space.
            % In the second example, a RIV stress, which requirs 3 inputs,
            % i.e. stage, conductance and river bottom elevation in that order
            % (see MODFLOW manual), the first spatial array DEM will be
            % interpreted as the spatial elevation relative to which the
            % values int the spreadsheet column 'H|stage' will be
            % taken. The second arguments, i.e., [], is empty and is
            % interpeted as 1, the scalar value 1 could also be supplied.
            % The tird argument, again DEM, will be taken as the plane of
            % elevations relative to which the column 'Hb*' in the
            % spreadsheet will be interpreted.
            % Leaving out arguments is the same as specifying []. To
            % specify only the last DEM the call should read:
            %
            % area2 = areaObj(gr,basename,sheetNm,'RIV',  [],[],DEM)
            %
            % alternatively one can explicitly associate an argument with a
            % column in the spreadsheet as follows
            %
            % riv = pointObj(gr,basename,sheetName,g'RIV',{DEM 'stage'},[],{DEM hBot})
            % or
            %
            % riv = pointObj(gr,basename,sheetName,g'RIV',{DEM 'stage' 'hBot'})
            %
            % in the latter case the columns 'stage' and 'hBot' must be
            % present int he worksheet that is read.
            %
            % EXAMPLES:
            %       riv = pointObj(gr,basename,sheetName,'RIV',{DEM 'stage' hBot},c)
            %       drn = pointObj(gr,basename,sheetName,'DRN',DEM,c})
            %       ghb = pointObj(gr,basename,sheetName,'GHB',DEM,c)
            %       chd = pointObj(gr,basename,sheetName,'CHD',{DEM 'Hstart' 'Hend'})
            %       flux= pointObj(gr,basename,sheetName,'WEL',q)
            %
            % Overview:
            %  STRESS         requiredHeadrs(1, 2 or 3)
            %                      1      2      3
            %  WEL|FLUX         flux|q
            %   DRN             H|elev*   C
            %   GHB             H|stage   C
            %   RIV             H|stage   C    Hb*
            %   CHD             H1|Hstart H2|Hend            
            %
            % instead of baseName,sheetName one may use 'struct',struct,
            % such as in:
            %       riv = pointObj(gr,'','','RIV','struct',struct,...)
            % That is, the fields basename and sheetName must be empty
            % strings and the parameter pair 'struct',struct must be
            % present in the argument list. In that case the fields in teh
            % struct will be used instead of the data info contained
            % normally in the sheetNm, in such a way, that the difference
            % will not be seen once the object is constructed. This method
            % allows generating objects from structs such as shapes,
            % avoiding the use of intermediate Excel tables. It is clear
            % that the fieldNames in the struct must be as required for
            % subsequent use of the objects.
            %
            %
            %
    %
    % TO 130818
    properties
        stressNames = {'NIL','WEL','FLUX','DRN','GHB','RIV','CHD'};
        name          % name of the object
        selected      % string used to select line from spreadsheet
        parent        % pointSeriesObj id to which point belongs
        type    % one of o.stressNames
        selectHdr     % name of header to based on which the lnes are selected
        activePerColName % name of column in PER sheet determinign when object is active (use 0 or 1)
        Idx           % global index of model grid cell(s) pertaining to this object
        A             % area of model grid cells pertainning to this object.
                      % For lineObj objects this is the intersection length
                      % between the line and the cell, not the area.        
        grSize        % size of the grid if it exists        
        vertexHdr     % headers of the inported vertex list (numerical values)
        vertex        % numeric vertex values as imported
        vertexTxtHdr  % headers of the imported vertex list (strings)
        vertexTxt     % char class vertex values imported
        iColx,iColy   % column numbers containing x and y coordinates
        narg          % number of arguments required by specific boundary type (see: modflow manual)
                      % i.e. WEL requires 1, RIV requires 3, other stresses require 2 
        V             % spatial values for each argument type V{1},V{2},V{3} obtained by
                      % interpolation of intersection center on point. There can be as many V
                      % as the value of narg (see narg)

        Dt; t;        % stress period length and time at end of stres period
        Q             % total flow into object
        CcolNames     % names of columns of species conc in order of species in model
                      % e.g.: {salinity temperature carbonate} or 'salinity'
        C             % input conc if object infiltrates
        Cout          % output conc of object. A weighted total of the mixed concentration
                      % of the species at end of each stress period. One
                      % species per line in order species are in input.
                      % Cout(3,4) would be Cout end of SP 4, species(3).
        NCOMP         % number of chemical species in simulation
                      % Cout=simulated dissolved concentration in point screen Q-mixed        

        wpix=2;       % pixelwidth of lines on screen for plotting        
        whdl=NaN(3,1);  % 3-valued handle: point casing and point screen in XSdrawing
        remark = '';  % any text
        
        FaceColor = 'w';  % outline color of plotted lines in XS or YS
        EdgeColor = 'b';  % screen coor of plotted lines in XS or YS
        FaceAlpha = 1;    % transparency of wel screen
        EdgeAlpha = 1;    % transparency of lines
        marker    = 'o';

        created      % moment at which point was create by the programm

        UserData  % all data not one of fixed fields and for any other purpose
    end
    properties (Dependent=true)
        x, y
    end
    methods
        function o=pointObj(varargin)
            %POINTOBJ: constructor of pointObj, see help pointObj
            % USAGE:
            %
            % obj =
            % pointObj([gr,]basename,sheetNm,stressType[selectionHdr,colName,][,opt,value,opt,value...,][DEM])
            %
            % basename of workbook with
            % sheet sheetname
            %   with the data in table form, i.e. a header line followed
            %   by data specifying the object(s).
            % stressType
            %   one of 'WEL','FLUX','DRN','GHB','RIV',CHD','NIL'
            % selectionColName,selectStr
            %   selectionColname is the name of column in worksheet to be used
            %   for deciding which lines to use. This is one by comparing the
            %   string values in the column with selectStr through function strmatchi.
            %   The obligatory 'name' column will subsequently be used to
            %   split the read lins into distinct objects.
            % opt.value,optvalue...
            %   pairs of option string with its value. Recognized pairs
            %   recognized are
            %   'active',perSheetColHdr
            %      uses the columns perSheetColHdr to decide during which
            %      stress periods this object is active (use 0 or 1 in that
            %      column to decideabout on or off.
            %      To use the opposit values in the same active column, place
            %      prepend the perSheetColHdrwith ~ (tilde). 
            %   'con',concentration
            %      to be used whith transport modelling
            % DEM is gr.Ny,gr.Nx plane of elevation data.
            %   DEM will be used to get elevation data from, unless that
            %   z-column is specified in the worksheet.
            
            if nargin<2, return; end
            
            [o.activePerColName,varargin] = getProp(varargin,'active','');
            [o.CcolNames       ,varargin] = getProp(varargin,'conc','');
            
            [myStruct, varargin] = getProp(varargin,'struct',{});
            
            if ~isempty(o.CcolNames) && ~iscell(o.CcolNames)
              o.CcolNames = {o.CcolNames};
            end          
          
            [gr,varargin] = getType(varargin,'gridObj',[]);
            
            isPointObj = strcmp(class(o),'pointObj');
            isAreaObj  = strcmp(class(o),'area2Obj');
            
            %% Get workbookName, sheetName and bcn/stress type            
            [basename    ,varargin] = getNext(varargin,'char','');
            [sheetNm     ,varargin] = getNext(varargin,'char','');
            [o.type      ,varargin] = getNext(varargin,'char','NIL');
            [o.selectHdr ,varargin] = getNext(varargin,'char','');
            
            o.type = upper(o.type);
            
            if isempty(basename) && isempty(myStruct)
                error('First argument must be the basename (type char, not <<%s>>.',...
                    class(varargin{1}));
            end
            if isempty(sheetNm) && isempty(myStruct)
                error('Second argument must be the sheetName');
            end
            if ~strmatchi(o.type,o.stressNames)
                error('Specified type <<%s>> not one of the stresses: \n<<%s>>',...
                    o.type,sprintf(' %s',o.stressNames{:}));
            end
            
            if ~isempty(o.selectHdr)
                [o.selected   ,varargin] = getNext(varargin,{'char','double','cell'},[]);
                if isempty(o.selected) ||...
                        (isnumeric(o.selected) && ~(isscalar(o.selected) || iscell(o.selected)))
                    error(['The fourth string <<%s>> should be the name of\n',...
                           'a column header in the worksheet and must be followed by\n',...
                           'another string indicating the value to select the lines\n',...
                           'from the data spreadsheet.\n'],o.selectHdr);
                end
            else
                o.selectHdr = 'name'; % default column to look
                o.selected      = '';     % will select all lines
            end

            %% Get the point vertices from myStruct or from basename/sheetname

            if ~isempty(myStruct)
                %% Get data from given struct and convert it so that
                %  the calling routine will not see the differece between
                %  calling with a struct or getting data from Excel
                myStruct =  getTableFromStruct(myStruct);
                
                ofields = fieldnames(o);
                sfields = fieldnames(myStruct);
                
                o = repmat(o,size(myStruct));

                mem = ismember(sfields,ofields);
                for io=numel(myStruct):-1:1
                    for i=numel(mem):-1:1
                        if mem(i)
                            field = sfields{i};
                            o(io).(field)=myStruct(io).(field);
                        end
                    end
                end
                
                [o.iColx]  = deal(strmatchi('x',o(1).vertexHdr,'exact'));
                [o.iColy]  = deal(strmatchi('y',o(2).vertexHdr,'exact'));
                [o.created]= deal(now);

                if isempty(gr),  return; end
                
                
            else
                %% Get data from Excel workbook
                
                [o.vertexHdr,o.vertex,o.vertexTxtHdr,o.vertexTxt] = getExcelData(basename,sheetNm,'Hor');

                % assert vertexHdr and vertexTxtHdr are unique
                % to prevent trouble with strmatchi not finding unique headers
                if numel(unique(o.vertexHdr)) ~= numel(o.vertexHdr)
                    error('Headers in vertexHdr are non unique in workbook <<%s>> sheet <<%s>>\n<<%s>>',...
                          basename,sheetNm,sprintf(' %s',o.vertexHdr{:}));
                end
                if numel(unique(o.vertexTxtHdr)) ~= numel(o.vertexTxtHdr)
                    error('Headers in vertexTxtHdr are non unique in workbook <<%s>> sheet <<%s>>\n<<%s>>',...
                        basename,sheetNm,sprintf(' %s',o.vertexTxtHdr{:}));
                end

                %% convert WGS coordinates to UTM if necessary and make sure the
                % x and y coordinatesa are the 2nd and 3rd column of vertex
                jlat = strmatchi('lat',o.vertexHdr);
                jlon = strmatchi('lon',o.vertexHdr);
                jx   = strmatchi('x'  ,o.vertexHdr,'exact');
                jy   = strmatchi('y'  ,o.vertexHdr,'exact');                        
                if (jx && jy)
                    % continune, skip
                else
                    if jx, o.vertex(:,jx) = []; end
                    if jy, o.vertex(:,jy) = []; end
                    if (jlat && jlon)
                            [xv yv] = wgs2utm(o.vertex(:,jlat),o.vertex(:,jlon));
                            o.vertex =    [ xv , yv ,o.vertex   ];
                            o.vertexHdr = ['x' ,'y',o.vertexHdr];
                            jx = 1;
                            jy = 2;
                    else
                        error('x,y or lat lon must exist as fields in sheet <<%s>>',sheetNm);
                    end
                end

                % handy for later use without searching
                o.iColx = jx;
                o.iColy = jy;

                %% select the lines to be used for the objects

                s = ['%s: there are no lines in column <<%d>> that\n',...
                      'correspond with the given property in excel file <<%s>> sheet <<%s>>'];

                % see if the the selectHdr is in vertexHdr or vertexTxtHdr
                jCol = strmatchi(o.selectHdr,o.vertexHdr);            
                if jCol % selection is based on numerical values (must be integers to work)
                    I = ismember(o.vertex(:,jCol),o.selected);
                    if ~any(I), error(s,mfilename,basename,sheetNm); end
                else % see if the selectHdr is in vertexTxtHdr                
                    jCol = strmatchi(o.selectHdr,o.vertexTxtHdr);
                    if jCol % selection is based on string
                        I = strmatchi(o.selected,o.vertexTxt(:,jCol));
                        if ~any(I), error(s,mfilename,jCol,basename,sheetNm); end
                    else
                        error('There is no column with header <<%s>> in excel file <<%s>> sheet <<%s>>',...
                            o.selectHdr,basename,sheetNm);
                    end
                end
                % select the lines
                if ~isempty(I)             
                    o.vertex    = o.vertex(   I,:);
                    o.vertexTxt = o.vertexTxt(I,:);
                end


                %% Split the input into as many lines as there are unique ids
                iColSelect = strmatchi(o.selectHdr,o.vertexHdr);                
                if iColSelect
                    isTxt = false;
                    Ids = unique(o.vertex(:,iColSelect),'first');
                else
                    isTxt = true;
                    iColSelect = strmatchi('name',o.vertexTxtHdr,'exact');
                    if ~iColSelect
                        iColSelect = strmatchi(o.selectHdr,o.vertexTxtHdr);
                    end
                    names = unique(o.vertexTxt(:,iColSelect),'first');
                end

                % remember b
                b   = o;         
                icx = o.iColx;
                icy = o.iColy;

                % split b into separate objects
                for io = numel(names):-1:1
                    o(io)          = b; 
                    if isTxt
                        o(io).name     = names{io};
                        I = ismember(b.vertexTxt(:,iColSelect),names{io});
                    else
                        o(io).name     = Ids(io);
                        I = ismember(b.vertext(:,iColSelect),Ids(io));
                    end
                    o(io).vertex   = b.vertex(   I,:);
                    o(io).vertexTxt= b.vertexTxt(I,:);
                    o(io).created  = now; % stamp

                    % CLOSE areaObj if required
                    if isAreaObj
                        allowableGap = 1;  % 1 m
                        xx = o(io).vertex([1 end],icx);
                        yy = o(io).vertex([1 end],icy);
                        if sqrt(diff(xx).^2+diff(yy).^2)>allowableGap
                            o(io).vertex(   end+1,:) = o(io).vertex(   1,:);
                           if ~isempty(o(io).vertexTxtHdr)
                                o(io).vertexTxt(end+1,:) = o(io).vertexTxt(1,:);
                           end
                        end
                    end
                end

                % Remove columns that are all NaN, allows using same table for
                % different puposes in Excel.
                for io=numel(o):-1:1
                    o(io).vertexHdr(:,all(isnan(o(io).vertex),1))=[];
                    o(io).vertex   (:,all(isnan(o(io).vertex),1))=[];
                end
            end  
            
            % Below we will have the grid available
            if isempty(gr),  return; end

            %% =======  Continue: Put Points into grid =====================

            % We are going to make sure that every object has either
            % relative z or abesolute z coordinates
            % on an object by object basis
            % To do this we need the grid


            for io = numel(o):-1:1

                o(io).grSize = gr.size;

                % Kick out columns with no data, i.e. NaN and abs(data)>1e10
                % This is not applicable to non-numeric columns
                J = all(abs(o(io).vertex)>=1e10,1);
                if any(J)
                    o(io).vertex(:,J)=[];
                    o(io).vertexHdr(J)=[];
                end


                %% Check to see that we have the correct info columns for this stress type

                allHdrs = [o(io).vertexHdr o(io).vertexTxtHdr];

                switch upper(o(io).type)
                    case 'NIL'
                        % skip
                        if ~strmatchi({'z','zRel'},o(io).vertexHdr)
                            o(io).vertex(:,end+1) = 0;
                            o(io).vertexHdr{end+1} = 'zRel';
                        end
                    case {'DRN','GHB'}
                        o(io).narg = 2;
                        if ~any(strmatchi( {'H','D'} ,allHdrs,'exact'))
                            error('Missing <<H&D>>  for %s in workbook <<%s>> sheet <<%s>>',o(io).type,basename,sheetNm);
                        end
                        isRel = strmatchi('D',allHdrs,'exact');
                    case 'RIV'
                        o(io).narg = 3;
                        if ~any(strmatchi( {'H','D'} ,allHdrs,'exact'))
                            error('Missing <<H&D>>  for %s in workbook <<%s>> sheet <<%s>>',o(io).type,basename,sheetNm);
                        end
                        if ~any(strmatchi( {'Hb','Db'}  ,allHdrs))
                            error('Missing bottom elev <<Hb&Db>> for %s in workbook <<%s>> sheet <<%s>>',o(io).type,basename,sheetNm);
                        end
                        isRel = strmatchi('Db',allHdrs,'exact');
                    case 'CHD'
                        o(io).narg = 2;
                        if ~any(strmatchi( {'Hst','H1','Dst','D1','H','D'} ,allHdrs))
                            error('Missing <<Hst&H1&D1&Dst&H|D>> for %s in workbook <<%s>> sheet <<%s>>',o(io).type,basename,sheetNm);
                        end
                        if ~any(strmatchi( {'Hen','H2','Den','D2','H','D'},allHdrs,'exact'))
                            error('Missing <<Hend&H2&D2&Dend&H&D>> for %s in workbook <<%s>> sheet <<%s>>',o(io).type,basename,sheetNm);
                        end
                end

                switch o(io).type
                    case {'WEL','FLUX'}
                        o(io).narg = 1;
                        if ~any(strmatchi({'q','flux'},allHdrs))
                            error('Missing data for ''q'' or ''flux''');
                        end                        
                    case {'DRN','GHB','RIV'}
                        if ~any(strmatchi( 'C'   ,allHdrs,'exact'))
                            error('Missing spec. conductance <<C>> for %s in workbook <<%s>> sheet <<%s>>',o(io).type,basename,sheetNm);
                        end
                end

                % The cell values to be written to the Modflow stress files
                % Will be separately and explicitly save in the object in
                % the arrays o(io).V which can have 1, 2 or 3 columns,
                % depending on the stress type.
                % The constructor accepts extra inputs that have the shape
                % gr.Ny,gr.Nx which will be used to modify the input from
                % the spredsheet. The meaning of the planes depends on the
                % stress type. If the respective stresstype data arbument is an
                % elevatoin, the data will be added to the plane (i.e. head
                % elevaton and river bottom). If not it will
                % be multiplie by the plane (i.g. the flux or conductance).
                % I expect this option ot te be used often but it can be
                % handy if for instance the conductance has to be
                % multiplied by a spatial multiplyer array.
                switch upper(o(io).type)
                    case {'WEL','FLUX'}, o(io).V = {1};
                    case {'DRN','GHB'},  o(io).V = {0 1};
                    case 'RIV',          o(io).V = {0 1 0};
                    case 'CHD',          o(io).V = {0 0};
                end

                if isPointObj

                    % Remove the points outside the model
                    msgId = 'list2short:ObjRemoved';
                    warning('on',msgId);

                    % Assert all vertices are inside the model boundaries
                    xv = o(io).vertex(:,o(io).iColx);
                    yv = o(io).vertex(:,o(io).iColy);
                    IDx = xyzindex(xv,yv,gr);  % I d'ont have Z avaialable here.
                                                    % indices in plane only !!
                    % remove vertices that are out size the x and y boundaries
                    % of the model
                    I = ~isnan(IDx);
                    if isempty(I)
                        warning(msgId,'%s %s removed, out of model!',class(o(io)),o(io).name);
                        o(io)=[];
                        continue
                    else
                        o(io).vertex    = o(io).vertex(I,:);
                        o(io).vertexTxt = o(io).vertexTxt(I,:);
                        xv  = xv(I);
                        yv  = yv(I);
                    end

                    % Check to see z is inside the model and which z-column to use
                    % So that we can compute the o(io).Idx

                    jZ  = strmatchi('z'   ,o(io).vertexHdr,'exact');
                    jrZ = strmatchi('zRel',o(io).vertexHdr');

                    needZ = strmatchi(o(io).type,{'WEL','FLUX','CHD','GHB'});

                    if needZ &&  ~(jZ|| jrZ)
                            error('%s %s requires z or zRel column',class(o(io)),o(io).name);
                    end

                    if jZ  % if we have a z-column
                        o(io).Idx = xyzindex(xv,yv,o(io).vertex(:,jZ),gr);

                    elseif jrZ % if we have a zRel column
                        % points above the model will be placed in layer 1
                        % points below the model will be placed in the
                        % bottom layer
                        iL = min(gr.Nlay,max(1,ceil(o(io).vertex(:,jrZ))));
                        o(io).Idx = xyzindex(xv,yv,gr) + (iL-1)*gr.Nxy;

                    elseif strcmpi(o(io).type,'DRN');

                        % if neither columns, be we have a DRN type allowing H
                        % or Hb or D or Db to be interpreted as z
                        if isRel
                            zv = interp2(gr.Xc,gr.Yc,gr.Z(:,:,1), xv, yv) - ...
                                o(io).vertex(:,strmatchi('D',o(io).vertexHdr,'exact'));
                            o(io).Idx = xyzindex(xv,yv,zv,gr);
                        else
                            zv = o(io).vertex(:,strmatchi('H',o(io).vertexHdr),'exact');
                            o(io).Idx = xyzindex(xv,yv,zv,gr);
                        end
                    elseif strcmpi(o(io).type,'RIV')

                        % same for RIV
                        if isRel
                            zv = interp2(gr.Xc,gr.Yc,gr.Z(:,:,1), xv,yv) - ...
                                o(io).vertex(:,strmatchi('Db',o(io).vertexHdr));
                            o(io).Idx = xyzindex(xv,yv,zv,gr);
                        else
                            zv = o(io).vertex(:,strmatchi('Hb',o(io).vertexHdr));
                            o(io).Idx = xyzindex(xv,yv,zv,gr);
                        end
                    end

                    % remove points outside the model (above or below)
                    I = find(~isnan(  o(io).Idx));
                    o(io).vertex    = o(io).vertex(I,:);
                    o(io).vertexTxt = o(io).vertexTxt(I,:);
                    o(io).Idx       = o(io).Idx(I);

                    % only for pointObj. Dont throw away vertices for other
                    % object types.
                    if isempty(I)
                        warning(msgId,' %s %s is empty, it is removed',class(o(io)),o(io).name);
                        o(io)=[];                        
                        continue;
                    end

                    % what if no objects are left ??
                    if isempty(o)
                        warning(msgId,'No objects of class %s left, probably all vertices are outside the model',class(o(io)));
                        continue
                    end

                    warning('off',msgId);


                    o(io).A = 1;  % gr.AREA3(o(io).Idx);

                    % Processing further input arguments.
                    % These are assumed to be planes of size [Ny,Nx]

                    % Set the planes and select the cells from them
                    % dependent on whether we have a pointObj or an area2Obj

                    for iPlane = 1:min(o(io).narg,numel(varargin))
                        if ~isempty(varargin{iPlane})

                            % If specified it must be a plane of the size of a layer
                            if ~all(size(varargin{iPlane}(:,:,1))==[gr.Ny,gr.Nx])
                                error('%s: argument must be of size Ny,Nx = [%d,%d]',mfilename,gr.Ny,gr.Nx);
                            end

                            o(io).V{iPlane} = varargin{iPlane}(o(io).Idx);

                            % Guarantee column vector
                            o(io).V{iPlane} = o(io).V{iPlane}(:);                    
                        end
                    end
                end
            end
        end

        function Cout=get.Cout(o), if isempty(o.Cout), Cout=NaN(size(o.C)); else Cout=o.Cout; end; end
        function C   =get.C(o),    if isempty(o.C),    C   =NaN(size(o.Q)); else C   =o.C;    end; end

        function o = check(o,IBOUND)
            %CHECK: remove cells that conincide with IBOUND==0
            mesId = [class(o) ':pointsAtZeroIBOUND'];
            warning('on',mesId);
            for io=numel(o):-1:1
                if nargin<2 || ~all(size(IBOUND)==o(io).grSize)
                    error('%s: need IBOUND in call and its size must be %s',...
                        mfilename,sprintf(' %d',o(iL).size));
                end
                nEl1 = numel(o(io).Idx);
                J = find(IBOUND(o(io).Idx)~=0);
                o(io).vertex    = o(io).vertex(J,:);
                o(io).vertexTxt = o(io).vertexTxt(J,:);
                for i=1:numel(o(io).V)
                    if numel(o(io).V{i})>1
                        o(io).V{i} = o(io).V{i}(J);
                    end
                end
                o(io).A         = o(io).A(J);
                o(io).Idx       = o(io).Idx(J);
                nEl2 = numel(o(io).Idx);
                if nEl1~=nEl2
                    warning(mesId,'%s: %s type = %s, name = %s, numel cells %d > %d (%d removed)',...
                        mfilename,class(o),o(io).type,o(io).name,nEl1,nEl2,nEl1-nEl2);
                end
                if isempty(o(io).Idx)
                    warning(mesId,'%s: %s type = %s, name = %s is empty',...
                              mfilename,class(o),o(io).type,o(io).name);
                    warning(mesId,'off');
                    o(io) = [];
                end
            end
            warning('off',mesId);
        end

        function plot(o,varargin)
            %PLOT: plots the vertices of the object(s) in 2D
            %
            % EXAMPLE: are2Obj.plot();
            %          pointObj.plot(ax);
            %          lineObj.plot(ax,'r','lineWidth',icx);
            %     ax is optional axis handle
            %     works for area2Obj, lineObj and pointObj
            %
            % TO 130825
            
            [ax,varargin] = getType(varargin,'axis',gca);
            
            isLspec = ~isempty(varargin) && isLineSpec(varargin{1});
            if ~isLspec, varargin = ['ro',varargin]; end
            
            for io=1:numel(o)
                if ~isLspec
                    varargin{1} = [mf_color(io), 'o'];
                end
                plot(ax,o(io).vertex(:,o(io).iColx),o(io).vertex(:,o(io).iColy),varargin{:});
            end
        end
        
        function plotCells(o,varargin)
            %PLOTCELLS: plot the cells pertaining to this object
            %
            % EXAMPLE:  area2Obj.plotCells();
            %           pointObj.plotCells(ax,'r');
            %           lineObj.plotCells('oc-','markerSize',5);
            %     ax is optional axis handle
            %     works for area2Obj, lineObj and pointObj
            %            
            % TO 130826
            
            for io=1:numel(o)
                if ~isempty(varargin) && isLineSpec(varargin{1})
                    lineSpec = varargin{1};
                else
                    if ~isempty(varargin), varargin(1)=[]; end
                    lineSpec = [o(io).marker o(io).EdgeColor o(io).marker];
                end
                plot3(gr.XM(o(io).Idx),gr.YM(o(io).Idx),gr.ZM(o(io).Idx),lineSpec,varargin);
            end
        end
        
        function fill(o,varargin)
            %FILL: plot H of pointObj plot as a circular marker with color
            % according to head field as shown in colorbar
            %
            % EXAMPLE: pointObj.fill();          % marker size will be 5
            %          pointObj.fill(markersize) % different marker size
            %          pointObj.fill(ax,markersize)
            %          pointObj(4:8).fill(ax,markersize)        % objects 4:8 only
            %
            % TO 130825
            
            defaultMarkerSize = 10;
            
            [ax , varargin] = getNext(varargin,'axis',gca);
            [markerSize, ~] = getNext(varargin,'double',defaultMarkerSize);
            
            for io=1:numel(o)
                % value to be plotted
                jCol = strmatchi('H',o(io).vertexHdr);
                h    = o(io).vertex(1,jCol);
                
                % get color from colormap
                clim = get(gca,'clim');
                cmap = colormap;
                iclr = max([1 ceil((h-clim(1))/diff(clim)*size(cmap,1))]);
                iclr = max(0, min(size(cmap,1), iclr));
                % plot circular marker with this color
                plot(ax,o(io).vertex(:,o(io).iColx),o(io).vertex(:,o(io).iColy),...
                    'o','markerSize',markerSize,'markerFaceColor',cmap(iclr,:));
            end
        end
        
        function o = label(o,varargin)
            %LABEL: plots label next to object at frequency of vertex
            %
            % EXAMPLE: pointObj.label('text');
            %          lineObj.label(ax,'lineLabel');
            %          areaObj.label(ax,[1 7 10],'areaLabel','color','r','fontsize',15,'angle',30);
            %
            % TO 130826
            
            [ax,varargin] = getNext(varargin,'axis',gca);
            [I,varargin]  = getNext(varargin,'double',1:10:1000);
            
            for io=1:length(o)
                idx = I(I<size(o(io).vertex,1));
                for i=idx(:)'
                    text(o(io).vertex(i,o(io).iColx),o(io).vertex(i,o(io).iColy),['\leftarrow ' sprintf('%s [%s]',o(io).name,o(io).type)],...
                        'verticalAlignment','baseline','rotation',30,'parent',ax,varargin{:});
                end
            end
        end        

        function o = label3(o,varargin)
            %LABEL: plots label next to object at frequency of vertex
            %
            % EXAMPLE: pointObj.label('text');
            %          lineObj.label(ax,'lineLabel');
            %          areaObj.label(ax,[1 7 10],'areaLabel','color','r','fontsize',15,'angle',30);
            %
            % TO 130826
            
            [ax,varargin] = getNext(varargin,'axis',gca);
            [I,varargin]  = getNext(varargin,'double',1:10:1000);
            
            for io=1:length(o)
                idx = I(I<size(o(io).vertex,1));
                for i=idx(:)'
                    text(ax,o(io).vertex(i,o(io).iColx),o(io).vertex(i,o(io).iColy),sprintf('% s [%s]',o(io).name,o(io).type),varargin{:});
                end
            end
        end
          
        function [Q,o] = Qtype(o,varargin)
            %QTYPE computes the injection [L2/t] of the object using the
            % the budget struct for the cells covered or intersected by
            % the object.
            %
            % USAGE:
            %   [q,lineObj] = lineObj.Qtype(B[,type][,'spec'])
            %   B is budget struct read by readBud(), B(end) is used
            %   axis is a valid axis
            %   type can be WEL,DRN,RIV,GHB,CHD
            %   if omitted, the type of the lineObj will be used.
            %   if axis given, the flows will be plotted as circles on that
            %   axis.
            %   q is a struct with the computed data for all objects.
            %   second output argument is lineObj with q stord in its
            %   UserData.
            %   'spec' signals to compute specific discharge of object
            %   which for pintObj = m3/d, lineObj m2/d and aea2Obj=m/d.
            %
            %   SEE ALSO  obj.setQ (more or less equivalent)
            %
            % TO 131012
            
                       % Budget labels corresponding to o(io).stressNames
            labels =      {'' 'WELLS' 'WELLS' 'DRAINS' 'HEADDEP'  'RIV'  'CONSTANTHEAD'};


            [spec,varargin] = getWord(varargin,'spec');
            [type_,varargin] = getType(varargin,'char','');
            [B ,         ~] = getType(varargin,'struct',[]);
            
            noType = isempty(type_);
            
            % Make sure varargin as at least a dummy lineSpec
            Nt = numel(B);
            
            for io=numel(o):-1:1                
                if noType, type_ = o(io).type; end                
                iStr = find(ismember(o(io).stressNames,type_) ,1,'first');
                label= labels{iStr};
                fldNm = ['Q' type_];
                o(io).UserData.(fldNm) = zeros( numel(o(io).Idx), Nt);
                for it=1:Nt
                    iTerm = find(ismember(labels,label),1,'first');
                    if iTerm
                        o(io).UserData.(fldNm)(1,it) = sum(B(it).term{iTerm}(o(io).Idx));
                    end
                end
                if isfield(B,'time')
                    o(io).UserData.(['t' type_]) = [B.time];
                end
                    Q = o(io).UserData.(fldNm);
                if spec
                    Q = Q./sum(o(io).A(:));
                end
            end            
        end  

        function o = setQ(o,varargin)
           %SETQ: computes and stores total disharge of this object (infiltration
           %positive as is always the case in MODFLOW
           %
           % USAGE: lineObj = lineObj.setQ(B[,type])
           %   B is the structure as obtained by readBud
           %   the output lineObj is required to store Q in it.
           %
           % SEE ALSO lineObj.printQ
           %
           % TO 131012
           
           [B,varargin ] = getType(varargin,'struct',[]);
           [type_,~    ] = getType(varargin,'char'  ,o(1).type); 
           
           if isempty(B), error('%s, needs Q, so call areaObj.printQ(B)',mfilename); end
          
           if nargout==0, error('linObj output argument is required to set Q'); end
           
           % Budget labels corresponding to o(io).stressNames
           labels =      {'' 'WELLS' 'WELLS' 'DRAINS' 'HEADDEP'  'RIVER'  'CONSTANTHEAD'};

           for io = numel(o):-1:1
                            
              if isempty(type_)
                  type_ = upper(o(io).type);
              else
                  type_ = upper(type_);
              end
              
              fldNm = ['Q' type_];

              iStr = find(ismember(o(io).stressNames,type_),1,'first');
              
              if isempty(iStr), error('stress type <<%s>> not recognized',type_); end

              if iStr==1
                  warning(msgId,'%s %s [%s], no discharge computed',class(o(io)),o(io).name,type_);
                  continue;
              end
              
              label = labels{iStr};              
              
              o(io).UserData.(fldNm) = zeros(numel(o(io).Idx),numel(B));

              for it = numel(B):-1:1
                  iTerm= strmatchi(label,B(it).label);
                  if iTerm
                      o(io).UserData.(fldNm)(:,it) = B(it).term{iTerm}(o(io).Idx);
                  else
                      % no cells of this stress in this stress period
                      o(io).UserData.(fldNm)(:,it) = 0;
                  end
              end
              
              if isfield(B,'time')
                  o(io).UserData.(['t' type_]) = [B.time];
              end

           end
        end
        
        function o = setH(o,varargin)
            %POINTOBJ/GETH -- retrieves the heads from the computed heads file
            % and ads time and hd field to UserData
            %
            % USAGE: pointObj = pointObj.setH([gr,]H[,userWantsRel])
            %
            %  gr is gridObj, if omitted gets only cell heads of pointObj
            %       if present, heads are exactly interpolated within their
            %       layer.
            %  H is heads structure as read by readDat: H(it).values(Ny,Nx,Nz)
            %  userWantsRel is either true or false telling whether absolute heads or
            %  are desired or relative with respect to ground surface elevation.
            %  Default for isRel is false.
            %
            % TO 131230

            if ~strcmpi(class(o),'pointObj')
                error('Works only for pointObj class');
            end
            
            [gr, varargin] = getType(varargin,'gridObj',[]);
            [H ,varargin ] = getType(varargin,'struct',[]);
            
            if isempty(H), error('Heads struct (as read by reaDat(...) is required'); end
            
            % isRel true means the user wants the head to be plotted
            % relative to DEM and not in abesolute elevations.
            [userWantsRel,~]= getNext(varargin,{'logical','double'},false);
            if userWantsRel && isempty(gr)
                error('grid input is required if groundwater depth is desired instead of head');
            end
            
            time = [H.time];
            Nt = numel(time);
            hd =NaN(numel(o),Nt);
            
            for io=numel(o):-1:1
                Nxy       = o(io).grSize(1)*o(io).grSize(2);
                xP(io)    = o(io).vertex(1,o(io).iColx);
                yP(io)    = o(io).vertex(1,o(io).iColy);            
            end
            
            IlayObj  = floor(([o.Idx]-1)/Nxy)+1;
                
            %% Process all objects at once    
            for iL = 1:o(1).grSize(end)
                % There can only be one layer valid, but it may be
                % different one for each object. Therefore, we have
                % to loop over all layers and pick out the objects
                % that are in the current layer
                Io = IlayObj==iL;  % Io are the objects that are in the current layer
                if any(Io)
                    if isempty(gr)
                        % Use head from cell in which pointObj is
                        for it=Nt:-1:1
                            hd(Io,it) = H(it).values([o(Io).Idx]);
                        end
                    else
                        % Interp head within layer
                        try
                            for it=Nt:-1:1
                                    hd(Io,it) = interp2(gr.Xc,gr.Yc,H(it).values(:,:,iL),xP(Io),yP(Io));
                            end
                        catch ME
                            fprintf('%s\n',ME.message);
                            for it=Nt:-1:1
                                hd(Io,it) = interp1(gr.xc,H(it).values(1,:,iL),xP(Io));
                            end
                        end
                        if userWantsRel
                            dem = interp2(gr.Xc,gr.Yc,gr.Z(:,:,1),xP(Io),yP(Io));
                            hd = bsxfun(@minus,hd,dem(:));
                        end
                    end
                end
            end

            for io=numel(o):-1:1
                o(io).UserData.time = time;
                o(io).UserData.hd   = hd(io,:);
            end
        end

        function o = plotHead(o,varargin)
            %POINTOBJ/PLOTHEAD -- plots and gives the heads at pointObj
            %locations, considering them as piezometers
            %
            % USAGE: pointObj.plotHead([ax][,gr,]H[,userWantsRel][,'meas',plotOptions])
            %
            %  ax is axis, default is gca
            %  gr is gridObj, if omitted plots cell heads of pointObj
            %       if present, heads are exactly interpolated within their
            %       layer.
            %  H is heads structure as read by readDat: H(it).values(Ny,Nx,Nz)
            %  userWantsRel is either true or false telling whether absolute heads or
            %  are desired or relative with respect to ground surface elevation.
            %  Default for isRel is false, meaning absoute heads will be
            %  plotted.
            %  'meas' if specfied means that the measured heads will also
            %  be plotted. If omitted they are not plotted.
            %  plotOptions as accepted by MatLab's plot function.
            %
            % TO 131011

            if ~strcmpi(class(o),'pointObj')
                error('Works only for pointObj class');
            end
            
            [meas,varargin]= getWord(varargin,'meas');
            [ax, varargin] = getType(varargin,'axis',[]);
            [gr, varargin] = getType(varargin,'gridObj',[]);
            [H ,varargin ] = getType(varargin,'struct',[]);
            
            if isempty(H), error('head struct as read by readDat(...) is required'); end

            % isRel true means the user wants the head to be plotted
            % relative to DEM and not in abesolute elevations.
            [userWantsRel,varargin]= getNext(varargin,{'logical','double'},false);
                        
            if isempty(gr)
                o = o.setH(H,userWantsRel);
            else
                o = o.setH(gr,H,userWantsRel);
            end
            
            time = [H.time];
                                    
            isLsp = ~isempty(varargin) && isLineSpec(varargin{1});

            if isempty(ax)                
                figure('name','heads of point objects','pos',screenPos(0.75));
                ax = axes('nextplot','add','xGrid','on','yGrid','on');
                xlabel(ax,sprintf('t --> [tmin =%d, tmax = %d]',time([1 end])));
                ylabel(ax,'head [m]');
                if userWantsRel
                    title('head - gr.surf.elev. vs time for pointOjb');
                else
                    title('abs heads vs time for pointObj');
                end
            else
                ax=gca;
            end
            
            if ~isLsp
                if ~isempty(varargin)
                    varargin = ['dummy' varargin];
                else
                    varargin ={'dummy'};
                end
            end

            k=0; leg={};
            clrs = 'brgkmcy';
            for io=1:numel(o)

                iD = find(ismember(o(io).vertexHdr,'D'),1,'first'); if isempty(iD), iD=0; end
                iH = find(ismember(o(io).vertexHdr,'H'),1,'first'); if isempty(iH), iH=0; end

                if meas
                    if iH
                        if userWantsRel
                           h = o(io).vertex(1,iH) - gr.Z(gr.IdTop(o(io).Idx));
                        else
                           h = o(io).vertex(1,iH);
                        end
                    else
                        if iD
                           if isempty(gr)
                               error('Plotting relative piezometer depths requires grid');
                           end
                           if userWantsRel
                               h =                           -o(io).vertex(1,iD);
                           else
                               h = gr.Z(gr.IdTop(o(io).Idx)) - o(io).vertex(1,iD);
                           end
                        end
                    end
                end
            
                if ~isLsp
                    iL = floor(io/numel(clrs))+1;
                    varargin{1} = [mf_color(io,clrs) mf_linetype(iL)];
                end
                k=k+1;
                leg{k} = o(io).name; %#ok
                plot(ax,o(io).UserData.time,o(io).UserData.hd,varargin{:});

                if meas && (iD || iH)
                    k=k+1;
                    leg{k} = [o(io).name ', data avg']; %#ok
                    plot(ax,o(io).UserData.time,h * ones(size(o(io).UserData.time)),[mf_color(io),'.']);
                end
            end
            legend(ax,leg{:},1);            
        end
        
        function o = printQ(o,varargin)
            %PRINTQ: prints the total discharge of the objects (extr. negative)
            %
            % USAGE [lineObj =] lineObj.printQ([B][,stressType][,dimenion][,'spec'])
            %  B = struct as obtained by readBud
            %  StressType is oneof ['WEL', 'DRN', 'GHB', 'RIV', 'CHD']
            %  Dimension  is oneof [{'m3/d'}, 'm3/y', 'Mm3/a', 'l/s']
            %  'spec' signals that output should be per unit length along
            %  the object.
            %  If B is omitted prints out Q(t) contained in object is
            %    lineObj = lineObj.setQ(B,varargin{:})
            %  If lineObj is also output argument, then t and Q will be set
            %  in the object.
            %
            % SEE ALSO lineObj.setQ
            %
            % TO 130825
            
            if isempty(o), return; end
            
            [spec ,varargin ] = getWord(varargin,'spec');
            [B  ,varargin   ] = getType(varargin,'struct',{});
            [type_ ,varargin] = getType(varargin,'char'  ,''); 
            [dimension ,~   ] = getType(varargin,'char'  ,'');
            
            if ~isempty(B)
                o = o.setQ(B,varargin{:});
            end
            
            if isempty(dimension) || ismember(upper(dimension),o(1).stressNames)
                    % switch dimension and type_
                    dum = dimension;
                    dimension = type_;
                    type_ = dum;
            end
            if isempty(dimension)
                dimension = 'md';
            end
            
            if ~isempty(type_)
                iStr = find(ismember(o(1).stressNames,type_),1,'first');
                if isempty(iStr)
                    error('%s %s type %s not oneof { %s}',class(o(1)),o(1).name,type_,sprintf('%s ',o(1).stressNames{:}));
                end
            end
            typeEmpty = isempty(type_);

            % Handle dimension text depending on object type and request
            % for total flow [L3/T] or object class specific flow, i.e.
            % m3/d for pointObj, m2/d for lineObj and m/d for areaObj
            dims = {
                'm /d',1
                'm /h',1/24
                'M /y',365.24/1e6
                'M /a',365.24/1e6
                'm /y',365.24/1e6
                'm /a',325.24/1e6
                'l /s',1/86.4};
            
            k = find(ismember(dims(:,1),[dimension(1) ' /' dimension(end)]),1,'first');
            if isempty(k)
                error(['Unknown discharge dimension <<%s>>, use oneof\n',...
                            'm3/d m3/h Mm3/a Mm3/y m3/a m3/y l/s'],dimension);
            end

            f = dims{k,2};
            
            if k==size(dims,1)
                if spec
                    dimStr = 'l/s/m';
                else
                    dimStr = 'l/s';
                end
            else
                dimStr = dims{k};
                if spec
                    switch class(o)
                        case 'pointObj', dimStr(2)='3';
                        case 'lineObj' , dimStr(2)='2';
                        case 'area2Obj', dimStr(2)='1';
                    end
                else
                    dimStr(2)='3';
                end
            end
                        
            msgId =  'warning:printQ';
            warning('on',msgId);
            
            % Print it
            [o.t]    = deal(     [0 [B.time]]);
            [o.Dt]   = deal(diff([0 [B.time]]));

            fprintf('\n');
            fprintf('%32sAverage t = %g - %g d,',' ',o(1).t([1 end]));
            fprintf('    t [d]->'); fprintf(' %12g',o(1).t);
            fprintf('\n');
            
            for io=1:numel(o)
                
                if typeEmpty, type_ = o(io).type; end
                
                fldNm = ['Q' type_];

                if ~isfield(o(io).UserData,fldNm)
                    warning(msgId,'No fieldName <<%s>> in UserData of %s %s',type_,class(o(io)),o(io).name);
                    continue;
                end
                fprintf('%32s [%5s],',o(io).name, type_);

                fprintf(' %12g %4s',sum(f*sum(o(io).UserData.(fldNm),1).*o(io).Dt,2)/sum(o(io).Dt),dimStr);

                if spec
                    fprintf('| %4s: ',dimStr);
                    fprintf(' %12g',f*sum(o(io).UserData.(fldNm),1)./sum(o(io).A(:)));
                else
                    fprintf('| %4s: ',dimStr);
                    fprintf(' %12g',f*sum(o(io).UserData.(fldNm),1));
                end
                fprintf('\n');
            end
            warning('off',msgId);
        end
        
        function write(o,fid,iPer)
           %% point.write(fid,iPer) -- write this point to file for this stress period
           %  using in writeBCN.m called by mf_setup
           %  TO 120629
           for io=1:length(o)
               if ~isnan(o(io).Q(iPer)) && o(io).Q(iPer)~=0,
                   for iCell = 1:length(o(io).idx)
                       fprintf(fid,'%10d%10d%10d%10g\n',...
                           o(io).LRC(iCell,:),...
                           o(io).fQ(iCell)*o(io).Q(iPer).');
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
       
        function o = sort(o,colName)
            % pointObj = pointObj.sort(colName)
            % sorts the objects according to colName, which may refer to
            % a numeric column in vertex or an text column in vertexTxt
            % TO 131108
            
            ic = strmatchi(colName,o(1).vertexHdr);
            if ic
                for io = numel(o):-1:1
                    list(io) = o(io).vertex(1,ic);
                end
                [~,I] = sort(list);
                o = o(I);
            else
                ic = strmatchi(colName,o(1).vertexTxtHdr);
                if ic
                    for io = numel(o):-1:1
                        list{io} = o(io).vertexTxt{1,ic};
                    end
                    [~,I] = sort(list);
                    o = o(I);
                else
                    error('Can''t find column <<%s>> neither in ...\nheaders <<%s >>\nnor in ...\nheaders <<%s >>',...
                        colName,sprintf(' ''%s''',o(1).vertexHdr{:}),sprintf(' ''%s''',o(1).vertexTxtHdr{:}));
                end
            end
        end
       
    end
end

function [BCN,PNTSRC]=bcnLine(o,varargin)
    %GRIDOBJ/BCNLINE --- generates input for boundary conditions/stresses
    % for objects given by a list of line vertices (line). These
    % objects can represent rivers, drains etc.
    %
    % USAGE:
    %    BCN         =gr.bcnLine(basename,type,line,vals)
    %    [GHB,PNTSRC]=gr.bcnLine(basename,type,line,vals,conc,'SP',SP);
    %    see examples below
    %
    % OUTPUTS:
    %   BCN is a cell array, with one element per stress period. Each
    %     element contains a list with the values per cell generated
    %     from the input as required by the specific package in MODFLOW or MT3MDS such
    %     as RIV, DRN, GHB, CHD ...
    %   PNTSRC
    %     if nargout>1 PNTSRC is generated for use by MT3DMS or SEAWAT
    %     and has the same structure. conc are only needed when PNTSRC
    %     is requested.
    %
    % Each element of the BCN cell array consists of the list of
    % cells with values for the specific stress period:
    %    [stressPeriodNr L R C values] for BCN
    %    [stressPeriodNr L R C conc1 ITYPE conc1 conc2 etc] for PNTSRC.
    %
    % You can specify to use only a limited set of stress periods
    % by the argument SP containgin stress period numbers,
    % or the string value pair
    % 'SP',stressPeriodNrs
    %
    % INPUTS:
    %   basename is the name of the Excel workbook holding the PER worksheet
    %   type     is the type of boundary condition according to MODLFOW,
    %            i.e. oneof 'DRN','RIV','GHB','CHD'.
    %   Data     point data.
    %   line     are a list of vertex coordinates in the following possible formats
    %       [ x y z ]       % 3D line specified in 3D
    %      {[x y] hRel}     % 2D line in layer iLay, where iLay is the
    %      layer number used in all cells in the output BCN and PNTSRC
    %      iLay has nothing to do with the elevation of drains, river
    %      bottoms etc.
    %
    %   vals must be a cell array {vals1 vals2 vals3 ...}
    %      vals are values of the values to be passed to the modflow
    %      cells after having converted them from the vertex values specified.
    %      They are stress-type specific. For instance, DRN
    %      requires elevation and conductance, GHB requires head and
    %      conductance, RIV requires stage conductance and bottom
    %      elevation, CHD requires two heads etc. All elevations, heads
    %      and stages are considered to be valid for the respective
    %      vertices. All conductances are considered valid for the
    %      vertex as well but per m of river length. This allows
    %      computing the cell conductance by multiplying this value by
    %      the length of the river crossing that particular cell. This
    %      length is computed from the intersetion of the line with the
    %      grid.
    %      Each of vals1 vals2 vals3 ... can be a specified in the following
    %      ways
    %        a scalar
    %        a vector           % length must equal number of vertices
    %        a plane array      % size must equal that of a plane of model cells
    %        a string           % value must correspond to a header in PER sheet
    %        a cell {scalar string}
    %        a cell {vector string}
    %        a cell {plane  string}
    %      Interprestation of the non-cell input:
    %      a string always refers to the transient values in the column in
    %      the PER sheet whose header equals the string. This results in a
    %      dfiferent value for each stress period for all vertices and so for
    %      all cells that are crossed by the line.
    %      a scalar refers to a single value for all cells for all stress
    %      periods.
    %      a vector must be as long as the number of vertices of the line
    %      and results in a different value for each vertex but the same
    %      value for each stress periods. The vertex values will be
    %      interpolated to the cells that are intersected by the line.
    %      a plane array may be used as areal information. The values at
    %      the intersected cells will be taken from the corresponding cells
    %      in the plane.
    %
    %      Interprestation of the cell input: {value|vector|plane string}
    %      By specifying a value through a cell as in the last sentence, the
    %      user can specify input that varies both in space and time.
    %      The numerical value|vector|plane values will be added to the
    %      stransient values corresponding to the string in case the value
    %      refers to an elevation or head. The value will be multiplied
    %      with the transient values in case the value corresponds to other
    %      variables, essentially conductance or flow in case type='WEL'.
    %    EXAMPLE:
    %      WEL = gr.bcnLine(basename,'WEL',{[x,y],L},Q);
    %        WEL wil contain a cell array of length NPER (number of stress
    %        periods), each element of which specifies a list of WEL cells.
    %        Each line of this list contains the folowing items:
    %        [SPnumber L R C QL]
    %        Notice that the layer number L in the call results in L being
    %        the same value in the list. QL corresponds to the value Q
    %        multiplied by the length of the line through the corresponding
    %        cell. So Q must have dimension [L2/T] (river width times q in [L/T]).
    %      WEL = gr.bcnLine(basename,'WEL',{x,y],2},{plane,Q});
    %        same as above, but the value QA in the list is now the product
    %        of Q with the corresponding cell value in plane and the area
    %        of the cell. The layer is now 2.
    %      RIV = gr.bcnLine(basename,'RIV',{[x,y],1},{plane,'stage'},C,{plane -5});
    %        This specifies a river resulting in lists with the following
    %        cell information
    %        [SPnumber L R C stage LC hBot]
    %        The layer number L is 1 in this case. The stage is a head or
    %        elevation value. It equals the value corresponding to the
    %        stress period number in column 'stage' in the PER worksheet
    %        plus the value in plane that corresponds with the cell. By
    %        specifying plane to be the DEM and stage negative values
    %        results in the stage being below ground surface at
    %        corresponding depth.
    %        The value LC is the value of C multiplied by the length of the
    %        line through the corresponding cell. So C must have dimension
    %        m/d so that LC has m2/d.
    %        The value hBot is, again, an elevation. Thefore the input
    %        {plane -5} results in a value of 5 m below the elevation
    %        defined by the given plane.
    %   conc    are concentrations.
    %       They are only needed for transport modelling, i.e. when
    %       PNTSRC is requested as output.
    %       The concentration values obey the same format and rules
    %       as the stress values explained above.
    %   SP or 'SP',SP can be used to specify a subset of stress
    %   periods to be included. SP is a list of numbers larger than
    %   zero and smaller than or equal to the number stress
    %   periods.
    %
    % TO 120413
    % TO 130613 (wrapped around bcnZone)

    %% assert input types
    [basename,varargin] = getNext(varargin,'char',[]);
    if ~ischar(basename)
        error('%s: First argument must be the basename of the workbook',mfilename);
    end

    [type,varargin] = getNext(varargin,'char',[]);
    if ~ischar(type)
        error('%s: Second argument must be string indicating the stress type like ''CHD'',''DRN'' etc.',mfilename);
    end

    [line,varargin] = getNext(varargin,{'double','cell'},[]);
    if isnumeric(line)
        switch size(line,2)
            case 2, line = {line 1};
            otherwise
                line = line(:,1:3);
        end
    elseif iscell(line) && ~isnumeric(line{1}(1))
        error(['%s: Third argument must define a line\n',...
               'either a 3 column array [x,y,z,] or a  cell {[x,y] zRel}'],...
               mfilename);
    end

    [vals,varargin] = getNext(varargin,{'cell','double','char'},[]);
    if isempty(vals)
        error(['%s: Fourth argument must contain the type-specific values to be\n',...
            'preferably a cell containing the values\n',...
            'each of which can be a sclar, a vector with length equal to the\n',...
            'number of vertices of the line or a string or\n',...
            'a cell  containing a {scalar|vector string} where\n',...
            'the scalar|verctor is a multiplier to the transient values\n',...
            'obtained from the column in the PER sheet with header corresponding\n',...
            'to the string'],mfilename);
    elseif ~iscell(vals)
        vals = {vals};
    end

    if nargout>1
        [conc, varargin] = getNext(varargin,{'cell','double','char'},[]);
        if ~(iscell(conc) || isnumeric(conc))
            error(['%s: Fifth argument must conain the concentrations in the form of\n',...
                'a vector with the number of columns equal to the number of species\n',...
                'an array of size (number of vertices numberof species)\n',...
                'string(s) for each species referring to a column in\n',...
                'the PER worksheet\n',...
                'REMEDY: if arg is name of column in PER worksheet, make it a cell by putting it in { }'],mfilename);
        end
    end

    % get stress periods as an extra argument
    % This seems an undocumented feature:
    % If desired, one can specify use of only certain stress periods instead of all
    [SP,varargin] = getProp(varargin,'SP',[]);
    if isempty(SP)
        [SP, ~] = getNext(varargin,'double',[]);
    end

    % Intersect line with the grid, to get line objects that hold
    % the information of each line section within a single grid.
    
    if iscell(line)
        iLay = line{end};

        if ~isscalar(iLay) || iLay<1 || iLay>gr.Nlay
            error('%s: layer number must be a scalar between 1 and %d',mfilename,gr.Nlay);
        end
        line = line{1};
    
        if ~isvector(line(1,:)) || ~isnumeric(line)
            error('%s: line must be numerical and have at least 2 columns',...
                mfilename);
        end

        % Intersect line with grid
        P = o.lineObjects(line,hRel);
    else
        if ~isnumeric(line) || size(line,2)<3
            error('%s: line must have at least 3 columns',mfilename);
        end

        %% In case the line has more than 3 fields (x,y,h)
        % then fields 4:end are considered values to be transferred to Modflow

        if size(line,2)>3
            vals = [line(:,4:end) vals];
        end

        P = o.lineObjects(line);
    end
    P = P.mergeDoubles();

    % Generate a zone array with true or 1 to indicate the cells
    % that are intersected by the line
    zoneArray = o.const(0);
    zoneArray([P.idx]) = true;

    % Each value in vals can be a scalar, numerical vector, a
    % string or a cell { } containing {scalar|vector string} combination.
    % The rows in vals, are given and interpreted as follows:
    % if vals is scalar, all vertices get this value.
    % if vals is vector, it must equal the length of the line, and
    % its values are attributed to the vertices of the line.
    % if vals is a string, all vertices will have the same value at
    % any stress period, where the values per stress period is
    % obtained from the column in the PER worksheet that
    % corresponds with its header.
    % if vals is a cell that combines a {sclar string} or {vector
    % string} combination then all stress period values in the
    % column in sheet PER correspoinding with the string are
    % multiplied by the scalar or the vector. The vector length
    % must comply with the number of vertices of the line.
    %
    % concentrations obey the same rules.

    %% Prepair intepolation along the line
    % cumulative line length at vertices
    sL = [0; cumsum(sqrt(sum(diff(line,1,1).^2,2)))];

    %Cumulative length along the intersected cells
    sCell = [0 cumsum([P.L])];  

    % Cumulative length to cell centers
    sCell = 0.5*(sCell(1:end-1)+sCell(2:end))';

    if size(line,2)==2
        % if hRel (or zRel) was given, line has two columns, use hRel.
        fr = hRel - floor(hRel); % fraction, rel z-coordinate in layer
        elevation = (o.ZTlay([P.idx])-o.DZ([P.idx])*fr)';
    else % else interpolate z-values (3rd column of line)

        % Take h to be the values of the third column
        % If line was 2D and zRel was given
        h  = line(:,3);                            

        % Elevation at intersected cell centers by interpolation
        elevation = interp1(sL,h,sCell);
    end

    % Interpolate the values given for the line to the cell centers
    for iv=1:numel(vals)
        if isnumeric(vals{iv}) && numel(vals{iv})>1
            if numel(vals{iv})~= size(line,1)
                error('%s: number of vals{%d}=%d must be 1 or equal to the number of line vertices (%s)',...
                    mfilename,iv,numel(vals{iv}),size(line,1));
            else
                % interpolate vertex values to the cell centers.
                % sCell is cumulative length of line to cell centers.
                vals{iv} = interp1(sL,vals{iv},sCell);
            end
        elseif iscell(vals{iv})
            if ischar(vals{iv}{1})
                vals{iv} = fliplr(vals{iv});
            end
            if ~isnumeric(vals{iv}{1})
                error('%s: vals{%d} must be a numeric string combination',...
                    mfilename,iv);
            elseif numel(vals{iv}{1}) ~= size(line,1)
                error('%s: number of vals{%d}=%d must be 1 or equal to the number of line vertices (%s)',...
                    mfilename,iv,numel(vals{iv}{1}),size(line,1));
            else
                % interpolate vertex values to the cell centers.
                % sCell is cumulative length of line to cell centers.
                vals{iv}{1} = interp1(sL,vals{iv}{1},sCell);
            end
        end
    end

    % Where applicable, multiply the given valus that are per m length
    % by the length of the cell they intersect.
    try
    switch type
        case 'WEL', zoneVals = {true elevation vals{1}.*[P.L]'}; % vals==q*L
        case 'DRN', zoneVals = {true elevation vals{1}.*[P.L]'}; % vals==cDrn
        case 'GHB', zoneVals = {true elevation vals{1}.*[P.L]'}; % vals==cGHB
        case 'CHD', zoneVals = {true vals{1} vals{2}};
        case 'RIV', zoneVals = {true elevation vals{1}.*[P.L]' vals{2}};
            % vals = {cRIV, bottom elevation}
        otherwise
            fprintf('%s: unknown TYPE %d, action skipped\n',mfilename,type);
    end
    catch ME
        error('%s: %s\n,check that in type %s the first argument is numeric',...
            mfilename,ME.message,type);
    end

    if nargin<6
        if isempty(SP)
            BCN          = bcnZone(basename,type,zoneArray,zoneVals);
        else
            BCN          = bcnZone(basename,type,zoneArray,zoneVals,'SP',SP);
        end
    else
        if isempty(SP)
            [BCN,PNTSRC] = bcnZone(basename,type,zoneArray,zoneVals,conc);
        else
            [BCN,PNTSRC] = bcnZone(basename,type,zoneArray,zoneVals,conc,'SP',SP);
        end
    end            
end

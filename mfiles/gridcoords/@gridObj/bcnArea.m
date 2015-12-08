function [BCN,PNTSRC]=bcnArea(o,varargin)
    %GRIDOBJ/BCNAREA --- generates input for boundary conditions/stresses
    % for objects enclosed by a polyline given by vertices (line). These
    % objects can represent areas with rivers, drains etc.
    %
    % USAGE:
    %    BCN         =gr.bcnArea(basename,type,{line,hRel},vals[,conc[,'SP',SP]])
    %    GHB         =gr.bcnArea(basename,'GHB',{[x,y],hRel},{heads C},SP])
    %    [GHB,PNTSRC]=gr.bcnArea(basename,'GHB',{[x,y],hRel},{heads C},{Conc1 Conc2 Conc3}[,'SP',SP
    %       Conc1 Conc2 Conc3 conc for three species
    %    CHD         =gr.bcnArea(basename,'CHD',{[x,y],hRel},{head1 head2});
    %    RIV         =gr.bcnArea(basename,'RIV',{[x,y],hRel},{head  C  bottomElev},SP)
    %       head(*) and Conc(*) and C may also be a string referring to transient
    %       values in the PER worksheet with header corresponding to the
    %       string.
    %       Conductance C values are per m2, they are multiplied by the
    %       area of the cells inside the polygon.
    %       If C is the combination of a string and numerical values like
    %           {string array}
    %       then the transient values in the column in the PER sheet with
    %       header equal to the string are multiplied by the numerical
    %       values. These values can be a scalar, a vector equal the number
    %       of cells in the polygon or an array corresponding to the size
    %       of a layer in the model.
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
    %   line     are a list of vertex coordinates in the following possible formats
    %      {[x y] zRel}   % area in layer given by zRel (iLay = ceil(zRel))
    %
    %   zRel     is the relative z-coordinate. It starts with 0 at the
    %      top of the model and ends with Nz. Therefore, the zRel=1
    %      represents the bottom of layer 1 or the top of layer 2, which is the
    %      same elevation. zRel = 4.3 is 30% into layer 5 from its top.
    %
    %   vals must be a cell array {vals1 vals2 vals3 ...}
    %      vals are values of the values to be passed to the modflow
    %      cells. They are stress type specific. For instance, DRN
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
    %      Each of vals1 vals2 vals3 ... can be a scalar or a
    %      string or a cell like {scalar string}. For instance:
    %      { {4 'stage'} 2.3 'heads' }
    %      A scalar is the value for the vertex as is.
    %      A string refers to the header of a column in the PER
    %      sheet to get a value for each stress period.
    %      A scalar string combination in a separate cell, where the
    %      string refers to the header of a column in the PER sheet
    %      with transient data, and the scalar is a multiplyer of these
    %      data specific to the vertex.
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

    [line,varargin] = getNext(varargin,{'cell','double'},[]);
    if isempty(line)
        error('%s: Third argument must define a line {[x y], zRel}\n',...
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

    % Make sure the values are cells, to generalize the input
    if ~iscell(vals), vals = {vals}; end

    % Intersect line with the grid, to get line objects that hold
    % the information of each line section within a single grid.

    if ~iscell(line)
        line = {line,0.5};
    end
    
    % in this case line is assumed to be given as a cell array
    % with two elements, the first being a numerical array with
    % two columns, being the x,y vertices and the second being
    % the relative z-coordinate starting 0 at the top of the
    % model and end with Nz. 1 is the bottom of layer 1, 0.5
    % its center, 4.5 the enter of layer 5 (not 4).
    hRel = line{2}(1);
    if ~isvector(line{1}(1,:)) || ~isnumeric(line{1})
        error('%s: line must be numerical and have at least 2 columns',...
            mfilename);
    end
    line = line{1}(:,1:2);

    iLay = ceil(hRel);

    zoneArray = false(o.size);
    zoneArray(:,:,iLay) = inpolygon(o.Xm,o.Ym,line(:,1),line(:,2));
    
    Idx = find(zoneArray);
    
    if isempty(Idx)
        error('%s: no cells within polygon',mfilename);
    end
    
    %fr   = hRel - floor(hRel); % fraction, rel z-coordinate in layer
    %elevation = (o.ZTlay(Idx)-o.DZ(Idx)*fr)';

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

    % Compute values by multiplying string values by area
    for iv=1:numel(vals)
        if ischar(vals{iv})
            vals{iv} = [{1} vals{iv}];
        end
    end
    
    for iv=1:numel(vals)
        if isnumeric(vals{iv})
            if numel(vals{iv}) == 1
                vals{iv} = vals{iv} * ones(size(Idx));
            elseif sameSize(vals{iv},o.AREA)
                vals{iv} = bsxfun(@times,vals{iv},ones(1,1,Nz));
                vals{iv} = vals{iv}(Idx);
            elseif sameSize(vals{iv},zoneArray)
                vals{iv} = vals{iv}(Idx);
            elseif numel(vals{iv})==numel(Idx)
                continue
            else
                error(['%s: number of vals{%d}=%d must be 1,\n',...
                    'or equal to the number of cells in the polyline (%d)\n',...
                    'or equal to a sheet|layer of the model [%d %d]'],...
                    mfilename,iv,numel(vals{iv}),numel(Idx),size(o.AREA));
            end
            if iv==1 &&  strcmpi(type,'WEL') % then vals{1} == Q, not head
                % multiply conductance with cell area
                vals{iv} = vals{iv} .* o.AREA(Idx);
            end
            if iv==2 && ~strcmpi(type,'CHD') % then vals{2} == C not head
                % multiply conductance with cell area
                vals{iv} = vals{iv} .* o.AREA3(Idx);
            end
        elseif iscell(vals{iv})
            if ischar(vals{iv}{1})
                vals{iv} = fliplr(vals{iv});
            end
            if ~isnumeric(vals{iv}{1})
                error('%s: vals{%d} must be a numeric string combination',...
                    mfilename,iv);
            end
            if numel(vals{iv}{1}) == 1
                vals{iv}{1} = vals{iv}{1} * ones(size(Idx));
            elseif samesize(vals{iv}{1},o.AREA)
                vals{iv}{1} = bsxfun(@times,vals{iv}{1},ones(1,1,Nz));
                vals{iv}{1} = vals{iv}{1}(Idx);
            elseif sameSize(vals{iv}{1},zoneArray)
                vals{iv}{1} = vals{iv}{1}(Idx);
            elseif numel(vals{iv}{1})==numel(Idx)
                % do nothing
            else
                error(['%s: number of vals{%d}=%d must be 1,\n',...
                    'or equal to the number of cells in the polyline (%d)\n',...
                    'or equal to a sheet|layer of the model [%d %d]'],...
                    mfilename,iv,numel(vals{iv}),numel(Idx),size(o.AREA));
            end
            if iv==1 && ~strcmpi(type,'CHD')
                % multiply conductance with cell area
                vals{iv}{1} = vals{iv}{1}.*o.AREA3(Idx);
            end
        end
    end

    % At this point the vals correspond to to Idx
    
    % Where applicable, multiply the given valus that are per m length
    % by the length of the cell they intersect.
 
    % put zone number in front of vals
    vals = [{true} vals];
    
    if nargin<6
        if isempty(SP)
            BCN          = bcnZone(basename,type,zoneArray,vals);
        else
            BCN          = bcnZone(basename,type,zoneArray,vals,'SP',SP);
        end
    else
        if isempty(SP)
            [BCN,PNTSRC] = bcnZone(basename,type,zoneArray,vals,conc);
        else
            [BCN,PNTSRC] = bcnZone(basename,type,zoneArray,vals,conc,'SP',SP);
        end
    end            
end

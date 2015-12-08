classdef mpath_particleGroupObj
    %MPATH_PARTICLEGROUPOBJ
    % definition of automatiically generated starting locations for Modpath
    % we may use here exactly the call to mpath_startingLocObj by passing the
    % arguments
    properties
        name                 = 'groupNm?';
        grid                 = 1;  % Must always be 1 (option for future versions)
        inputStyle           = 3;  % Automatically defined
        releaseOption        = 3;  % Always use multiple release times
        releaseStartTime     = 0;
        releasePeriodLength  = 0;
        releaseEventCount    = 1;
        releaseTimes         = 0;   % Vector of releaseTimes for group
        placementOption      = 0;   % Automatically generated
        gridCellRegionOption = 0;   % Automatically generated
        CHeadOption          = 2;   % Particles are also generated in constant head cells
        maskLayer            = 0;   % Layer number if maskArray is one layer only
        mask                 = [];  % Mask array is a
                                    %    logical array indicating celles that belong to a group
                                    %    can be 3D (gr.size) or
                                    %    2D (gr.size(1:2)) in combination with maskLayer
                                    %    or can define gridBlocks[ix1 ix2 iy1 iy12 iz1 iz2]
                                    %    may hvae more than 1 layer to
                                    %    define multiple gridBlocks at once
        particles            = [];  % of type mpath_startingLocObj
        IFace                = 0;   % Vector of group release ifaces (either 0 or vector ifaces, all >0)
        placement            = [1 1 1]; % vector of subdivisions for placement in cell according to IFace
        lineSpec             = 'ro';
        endPoints            = []; % endpoints array
        endpColHdr           = []; % labels of columns of endpoint array
        tsrPoints            = []; % time series points array
        tsrColHdr            = []; % time series points hdr
        pathColHdr           = []; % pathl line column header
        pathPoints           = []; % path line points
    end
    methods
        function o=mpath_particleGroupObj(varargin)
            % INPUTSTYLES 1, 2 or 3
            %
            % InputStyle 1
            % The location, releasetime and particle ID are specified explicitly for each particle in the simulation.
            % This style provides the flexibility to generate starting locations based on particle output from a
            % previous MODPATH simulation and retain the same particle ID values from one simulation to the next.
            %
            % InputStyle 2
            % Spatial locations are specified explicitly for one or more particle groups. Particle release time information
            % is specified for each particle group and applied to all particles in a group. Both single and multiple
            % particle release times are supported. MODPATH generates the required particles at the specified locations
            % and release times, and particle ID values are generated automatically. Style 2 allows complete flexibility
            % in specifying the spatial location of particles but retain the convenience and compactness of specifying
            % release time characteristics by particle group. Input style 2 produces a smaller data file than input style 2.
            %
            % InputStyle 3
            % Particle starting locations are generated for user-specified regions of grid cells based on a specified
            % template of particle locations located on either cell faces or within the cell for all cells in the specified
            % regions. Both single and multiple particle release times are supported. MODPATH generates the required
            % particles at the specified locations and release times, and particle ID values are generated automatically.
            % Style 3 is equivalent to the format used to specify starting locations directly in the MODPATH simulation
            % file. Input style 3 produces a much smaller data file than either style 1 or style 2.
            %            
            % USAGE 1: general input structure (inputstyles 1,2, 3)          
            %   particleGroupObj(gr,zoneVals[,varNm,value,varNm,value,...])
            %   where zoneVals is
            %        {zones,varNm,value,varNm,value,...; ...
            %         zones,varNm,value,varNm,value,...; ...
            %         ...};
            %   in which each line corresponds to a particleGroupObj to be
            %   generated.
            %   The varNm,value pairs in the call are the default for all particleGrpObj to be generated.
            %   those on each line of zoneVals are specific to the
            %   particleGrpObj on which line they are specified.
            %   zones specifies the particles to be generated
            %     zones  logical 3D array defining cells to place particles
            %     zones  numeric array with 3 columns defining XYZ
            %     coordinates of particles
            %     zones is of clas mpath_particleObj defines particles
            %     zones is numeric array of width 6 defines 3D blocks of
            %     cells in the grid to define particles.
            %
            % USAGE
            %   pGrp = mpath_particleGroupObj(gr,{mpath_particleObj,varNm,value,varNm,value,...; ...},varNm,value,...)
            %   pGrp = mpath_particleGroupObj(gr,{XYZcoordinates   ,varNm,value,varNm,value,...; ...},varNm,value,...)
            %   pGrp = mpath_particleGroupObj(gr,{cellBlockDef     ,varNm,value,varNm,value,...; ...},varNm,value,...)
            %   pGrp = mpath_particleGroupObj(gr,{logicalZoneArray ,varNm,value,varNm,value,...; ...},varNm,value,...)
            %
            %   The first two usage lines define groups with specified particles based on coordinates
            %   The third lines species cells in terms of 3D bllocks [ix1 ix2 iy1 iy2 iz1 iz2],
            %   which corresponds to input style 3.
            %   The fourth line corresponds to input style 2 with a zone
            %   array specifying which cells pertain to a zone. This
            %   zoneArray must be a logical 3D array of size gr.size.
            %
            % EXAMPLES
            % USAGE:
            %    pgrpObj = mpath_particleGroupObj(gr,particles,'releaseTimes',[1 10 100],'name','duck lake');
            %    pgrpObj = mpath_particleGroupObj(gr,{particles,'releaseTimes',[1 10 100]},'name','duck lake');
            %    both are equivalent because only one particle group is
            %    specified.
            %
            %    pgrpObj = mpath_particleGroupObj(gr,{globalXYZ,'releaseTimes',[1 10 100],'name','duck lake'});
            %    with only one particleGroup the { }  may be omitted
            %
            %    particleGrps = {particles 'name','duck lake' ,...
            %                    globalXYZ,'name','goose lake'};
            %      coordinates may be mixed, as it is the same inputStyle
            %
            %    pgrpObj = mpath_particleGroupObj(gr,particleGrps,'releaseTimes',[1 10 100]);
            %    redefine releastimes for existing particle group objects
            %
            %    pgrpObj = mpath_particleGroupObj(gr,{zoneArray,'releaseTimes',0:100:1000});
            %      zoneArray must be of class logical and of size equal to gr.size
            %
            %    example of zoneVals for 3 zones
            %    zoneVals={zoneArray1,'name','duck lake' ,'IFace',[1 2 6],'placement',[3 3]);
            %              zoneArray2,  'name','goose home','IFace',0,      'placement',[4 4 4];
            %              zoneArray3,  'IFace',  [5 6],       'place',[2 5]}, 'group Name','cow meadows'};
            %
            %    pgrpObj = mpath_particleGroupObj(gr,gridBlocks,'releaseTimes',0:100:1000);
            %        where gridBlocks may be [ix1 ix2 iy2 iy2 iz1 iz2; ...]
            %
            % here you notice that
            %   1) each line in zonVals starts with a definition of the cells pertaining to the particleGroupObj
            %   2) each line in zoneVals has the same number of varNm,value pairs
            %   3) the order of the varNm,value pairs does not matter
            %   4) autocompletion works (see abbeviation of 'placement' to
            %      'place' in the 3rd line
            %
            % 'IFace' defines the ifaces on which the particles are to be
            % placed, while 'placement' defines how the particles are
            % distributed across the different ifaces. Always use 3 values
            % for IFace telling the number of particle columns, rows and
            % layers that particles will be placed.
            %
            % SEE ALSO
            %  mpath_startingLocObj
            %
            % TO 130202
            
            global basename
            
            if nargin==0, return; end
            
                        %% Get default properties from workBook
            if isempty(basename)
                error(['basename not visible in %s.\n',...
                    'To resolve this, add command: <<global basename>> without << >>\n',...
                    'before assignment of basename in mf_adapt or mf_build'],mfilename);
            end

            % Make sure the grid is present in the input
            [gr,varargin] = getType(varargin,'gridObj',[]);
            if isempty(gr)
                error('must provide gridObj');
            end

            % zoneVals (i.e. multiple group specification) should be given
            % using cell format. Without cell format, consider all
            % subsequent input as a zoneVals cell array for one group.
            [zoneVals ,varargin]   = getNext(varargin,'cell',[]);
            if isempty(zoneVals)
                zoneVals = varargin;
                varargin={};
            end
            
            % If varargin is still not empty, then it should contain a list
            % of varNm,value pairs that will set properties for all groups,
            % which may be overridden by the settings that are group
            % specific.
            vararginToDo = varargin;
            if ~isempty(vararginToDo)                            
                if ~all(cellfun(@ischar,vararginToDo(1:2:end),'UniformOutput',false));
                    display(vararginToDo)
                    error('inputs must be varNm,varValue pairs');
                else
                    % varargin={};
                end
            end

            %% Get the settings specific for this MPATH run from the MPATH worksheet
            [mpathNams,mpathVals] = getExcelData(basename,'MPATH','vertical');
            
            %% Defaults for mpath_particleGroupObj groups will be stored in a template
            template = mpath_particleGroupObj();
            template.releaseOption       = mpathVals(strmatchi('releaseOption',      mpathNams),1);
            template.releaseStartTime    = mpathVals(strmatchi('releaseStartTime',   mpathNams),1);
            template.releaseEventCount   = mpathVals(strmatchi('releaseEventCount',  mpathNams),1);
            switch template.releaseOption
                case 1
                case 2
                    template.releasePeriodLength = mpathVals(strmatchi('releasePeriodLength',mpathNams),1);
                case 3
                    template.releaseTimes = mpathVals(strmatchi('releaseTimes',mpathNams),:);
                    template.releaseTimes(isnan(template.releaseTimes))=[];
                    if numel(template.releaseTimes)>=template.releaseEventCount
                        template.releaseTimes=template.releaseTimes(1:template.releaseEventCount);
                    else
                        error(['%s: insufficient times <<%d>> in releaseTimes in worksheet,',...
                                'must be <<%s>>'],...
                                numel(template.releaseTimes),template.releaseEventCount);
                    end
                otherwise
                    error(['Unknown option <<%d>> for releaseOption,\n',...
                        'Use 1..3, see workshet MPATH'],template.releaseOption);
            end

            %%
            template.CHeadOption = mpathVals(strmatchi('CHeadOption',        mpathNams),1);            
            
            %% Generate the input for all groups
            % Notice that the input style is the same for all particle
            % groups. This is a limitation of MPATH.
            
            % We have to direct the course of this particle group
            % generation depending on the methods of input supported by
            % mfLab.
            
            NGrp = size(zoneVals,1);
            
            for iGrp = NGrp:-1:1
                                        % template
                o(iGrp) = template;
                % defaults for all groups

                o(iGrp) = setProps(o(iGrp),vararginToDo);

                if iscell(zoneVals{iGrp,1})
                    zoneDef = zoneVals{iGrp}{1};
                    zVals   = zoneVals{iGrp}(2:end);
                else
                    zoneDef = zoneVals{iGrp,1};
                    zVals   = zoneVals(iGrp,2:end);
                end
                
                switch class(zoneDef)
                    case 'mpath_particleObj'
                            o(iGrp).inputStyle = 1;
                            o(iGrp).particles  = zoneVals{iGrp,1};
                    case 'double'
                        switch size(zoneDef,2)
                            case 3 % zoneVals{iGrp,1} is globalXYZ
                                o(iGrp).inputStyle = 2;

                                % generate particles using other options in rest of varargin{1}{iGrp,2:end}
                                o(iGrp).particles = mpath_particleObj(gr,zVals,zVals{2:end});
                            case 6 % rectangular blocks of cells
                                o(iGrp).inputStyle = 3;
                                o(iGrp).gridCellRegionOption = 1;
                                o(iGrp).mask = zoneDef;
                            otherwise
                                error(['unknown option for zoneDef\n' ...
                                    'REMDY: zoneDef should be logical if full 3D zones are to be defined']);
                        end
                    case 'logical'  % using MASK
                        o(iGrp).mask = zoneDef;
                        if size(zoneDef,3)>1 && all(size(zoneDef)==gr.size)
                            o(iGrp).gridCellRegionOption = 2;
                        elseif all(size(zoneDef(:,:,1)==gr.size(1:2)))
                            o(iGrp).gridCellRegionOption = 3;
                            o(iGrp).maskLay = getProp(zVals       ,'maskLay',[]);
                            o(iGrp).maskLay = getProp(vararginToDo,'maskLay',o(iGrp).maskLay);
                            if isempty(o(iGrp).maskLay) && isdoulbe(zVals{1})
                                o(iGrp).maskLay = zVals{1};
                                zVals = zVals(2:end);
                            end
                        else
                            error('unknown option for zoneDef');
                        end
                end
                % overwrite general options with specific options
                o(iGrp) = setProps(o(iGrp),zVals{:});

                %% placementOption
                if numel(o(iGrp).placement)==1
                    o(iGrp).placement = o(iGrp).placement([1 1 1]);
                end

                %% Cleaning IFace and placement

                % remove possible duplicate IFace numbers
                o(iGrp).IFace = unique(o(iGrp).IFace);

                % make sure IFace==0 is not mixed with other IFace
                % values, because of the separate input required
                % for modflow6 (see manual of modflow6)
                if any(o(iGrp).IFace==0)
                    o(iGrp).placementOption = 2;
                    if numel(o(iGrp).IFace)>1
                        error(['%s: IFace==0 must be a separate mpath_particleGroupObj:\n',...
                            'Therefore, separate IFace==0 from other values of IFace'],mfilename);
                    elseif numel(o(iGrp).placement)<3
                        % guarantee placement is a 3 value vector
                        % in case IFace==0 (1 value would be ok but
                        % is dealt with above
                        error('%s: with IFace==0 placement must be a 1 or 3-value vector',mfilename);
                    end

                % make sure we have a correct 3-value vector for
                % each IFace
                else
                    o(iGrp).placementOption = 1;
                    if numel(o(iGrp).placement) < 3
                        if all(ismember(o(iGrp).IFace,[1 2]))
                            o(iGrp).placement = [o(iGrp).placement(1) o(iGrp).placement(2) 1];
                        elseif all(ismember(o(iGrp).IFace,[3 4]))
                            o(iGrp).placement = [o(iGrp).placement(1) 1 o(iGrp).placement(2)];
                        elseif all(ismember(o(iGrp).IFace,[5 6]))
                            o(iGrp).placement = [1 o(iGrp).placement(1) o(iGrp).placement(2)];
                        else
                            error(['%s: mpath_particleGroupObj(%d), name<<%s>>:\n',...
                                'For this mix of IFaces, use a 3-value placement vector [NL,NR,NC]\n',...
                                'to allow unique placement of particles on the specified IFaces.'],...
                            mfilename,iGrp,o(iGrp).name);
                        end
                    end
                end                        
            end
        end
        function  o = getParticles(o,gr)
            % p = particles(o) -- generate particles for particleGroupObj
                       
            for iGrp = 1:numel(o)
            
                switch o(iGrp).inputStyle
                    case {1,2}
                        % already has particles
                    case 3
                        o(iGrp).particles = ...
                            mpath_particleObj(gr,...
                            o(iGrp).mask,...
                            'gridCellRegionOption',o(iGrp).gridCellRegionOption,...
                            'IFace',o(iGrp).IFace,'placement',o(iGrp).placement,...
                            'name',o(iGrp).name,'groupNr',iGrp,...
                            'releaseTimes',o(iGrp).releaseTimes');
                    otherwise
                        error('%s: Unknown inputStyle <<%d>> for particleGroup <<%d>> name <<%s>>',...
                            mfilename,o(iGrp).inputStyle,iGrp,o(iGrp).name);
                end
            end
        end
        
        function plot(o,varargin)
            % particleGrp.plot(varargin) --- plot particles of particleGroup(s)
            % varargin are plot options like with the original plotting function.
            if ~isempty(varargin) && isLineSpec(varargin{1})
                for iGrp = 1:numel(o)
                    o(iGrp).lineSpec = varargin{1};
                end
                varargin(1)=[];
            end
            for iGrp=1:numel(o)
                XYZ = vertcat(o(iGrp).particles.globalXYZ);
                plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),o(iGrp).lineSpec,varargin{:})
            end
        end
        
        function writeStartingLocations(o,fnameOrFid)
            % mpath_particleGroup.write(fid,unit)
            % writes particle groups to a mpath input file (always inputStyle 2)
            
            ctrlRecordDesired = true; % for use in warray
            freeFormat        = true; % for use in warray
            
            writeToFile = ischar(fnameOrFid);
            
            if writeToFile
                fid=fopen(fnameOrFid,'wt');
                %0.
                fprintf(fid,'%s\n',['# MATLAB  ' mfilename ' ' datestr(now)]);
                fprintf(    '%s\n',['# MODPATH ' mfilename ' ' datestr(now)]);            
            else
                fid=fnameOrFid;
            end

            nGroup = numel(o);
                        
            INPUTSTYLE = o(1).inputStyle;
            for iGrp = 1:nGroup                
                if o(iGrp).inputStyle ~= o(1).inputStyle
                    error(['%s: inputStyle of group %d: <<%d>> differs from that of group 1: <<%d>>,\n',...
                        '  all groups must have the same input style.'],...
                        mfilename,iGrp,o(iGrp).inputStyle,o(1).inputStyle);
                end
            end

            %1
            if writeToFile
                fprintf(fid,'%d INPUTSTYLE\n',INPUTSTYLE);
            elseif INPUTSTYLE ~=3
                error(['%s: You set particleGenerationOption=1 in worksheet MPATH.\n',...
                    'This implies that you intend to specify the info for modpath6 to\n',...
                    'fully automatically generate particles right in the simulaton file and not\n',...
                    'in a separate particleStartingLocations file. However, the inputStyle\n',...
                    'that is implied from the info specified in the mpath_particleGroupObj is non\n',...
                    'automatic or semi-automatic as it uses eiher predefined particles or coordinates.\n',...
                    'To resolve this, set particleGenerationOption=2 in worksheet MPATH,\n',...
                    'so that the particle starting location (definitions) will be written to a separate file,\n',...
                    'i.e. the particleStartingLocationsFile.'],mfilename);
            else
                % skip as INPUTSTYLE==3
            end
            
            switch INPUTSTYLE
                case 1
                    ID = 1;
                    %2
                    fprintf(fid,'%d     groupCount\n',nGroup); % groupCount
                    for iGrp = 1:nGroup
                        %3
                        fprintf(fid,'%s     name of particle group %d\n',o(iGrp).name,iGrp);
                        %4
                        fprintf(fid,'%d     particleCount\n',numel(o(iGrp).particles)*numel(o(iGrp).releaseTimes));
                        %5
                        o(iGrp).particles.write(fid,ID,o(iGrp).releasTime,o(iGrp).name);
                        ID = ID+numel(o);
                    end
                case 2
                    %6
                    fprintf(fid,'%d     groupCount\n',nGroup); % groupCount
                    for iGrp = 1:nGroup
                        %7
                        fprintf(fid,'%s     name of particle group %d\n',o(iGrp).name,iGrp);
                        %8
                        fprintf(fid,'%d %g %d     locationCount releaseStartTime releaseOption\n',...
                            numel(o(iGrp).particles),...
                            o(iGrp).releaseTimes(1),...
                            o(iGrp).releaseOption);

                        o(iGrp).releaseEventCount = numel(o(iGrp).releaseTimes);
                        switch o(iGrp).releaseOption
                            case 1
                            case 2
                                %9
                                fprintf(fid,'%d %g     releaseEventCount releasePeriodLength\n',...
                                    o(iGp).releaseEventCount,...
                                    diff(o(iGrp).releaseTimes([1 end]))/(o(iGrp).releaseEventCount-1));
                            case 3
                                %10
                                fprintf(fid,'%d     releaseEventCount\n',o(iGrp).releaseEventCount);
                                %11
                                fprintf(fid,' %g',o(iGrp).releaseTimes(2:end));
                                fprintf(fid,'\n');
                            otherwise
                                error('%s: illegal releaseOption <<%d>> use 1..3',...
                                    mfilename,o(iGrp).releaseOption)
                        end
                        %12
                        o(iGrp).particles.writeLoc(fid);
                    end
                case 3
                    %13
                    fprintf(fid,'%d     groupCount\n',nGroup); % groupCount
                    for iGrp = 1:nGroup
                        %14
                        fprintf(fid,'%s     name of particle group # %d\n',o(iGrp).name,iGrp);
                        %15
                        fprintf(fid,'%d %d %d %g %d %d     grid gridCellRegionOption placementOption releaseStartTime releasOption CHeadOption\n',...
                            o(iGrp).grid,...
                            o(iGrp).gridCellRegionOption,...
                            o(iGrp).placementOption,...
                            o(iGrp).releaseTimes(1),...
                            o(iGrp).releaseOption,...
                            o(iGrp).CHeadOption);
                        o(iGrp).releaseEventCount = numel(o(iGrp).releaseTimes);
                        switch o(iGrp).releaseOption
                            case 1
                            case 2
                                %16
                                fprintf(fid,'%d %g     releaseEventCount releasePeriodLength\n',...
                                    o.releaseEventCount,...
                                    diff(o(iGrp).releaseTimes([1 end]))/(o(iGrp).releaseEventCount-1));
                            case 3
                                %17
                                fprintf(fid,'%d     releaseEventCount\n',o(iGrp).releaseEventCount);
                                %18
                                fprintf(fid,' %g',o(iGrp).releaseTimes);
                                fprintf(fid,'\n');
                            otherwise
                                error('%s: illegal releaseOption <<%d>> use 1..3',...
                                    mfilename,o(iGrp).releaseOption)
                        end
                        switch o(iGrp).gridCellRegionOption
                            case 1
                                %19
                                % in this case o(iGrp).zone = [minLay minRow minCol maxLay maxRow maxCol]
                                fprintf(fid,'%d %d %d %d %d %d     minLay minRow minCol maxLay maxRow maxCol\n',o(iGrp).mask);
                            case 2
                                %20
                                for j=1:size(o(iGrp).mask,3)
                                    warray(fid,o(iGrp).mask(:,:,j),0,'(40I2)','Mask(NCOL,NROW)',ctrlRecordDesired,freeFormat);
                                end
                            case 3
                                %21
                                mLay = o(iGrp).maskLayer;
                                fprintf(fid,'%d     maskLayer\n',mLay);
                                %22
                                warray(fid,o(iGrp).mask(:,:,1),0,'(40I2)','Mask(NCOL,NROW)',ctrlRecordDesired,freeFormat);
                            otherwise
                                error('%s: illegal gridCellRegionOption <<%d>> use 1..3',...
                                    mfilename,o(iGrp).gridCellRegionOption);
                        end
                                                
                        switch o(iGrp).placementOption
                            case 1
                                %19
                                fprintf(fid,'%d     faceCount\n',numel(o(iGrp).IFace));
                                for iface = o(iGrp).IFace(:)'
                                    %20
                                    switch iface
                                        case {1,2}
                                            fprintf(fid,'%d  %d %d     IFace  NLay NRow\n',iface,o(iGrp).placement([1 2]));
                                        case {3,4}
                                            fprintf(fid,'%d  %d %d     IFace  NLay NCol\n',iface,o(iGrp).placement([1 3]));
                                        case {5,6}
                                            fprintf(fid,'%d  %d %d     IFace  NRow NCol\n',iface,o(iGrp).placement([2 3]));
                                        otherwise
                                            error('%s: illegal iface %d',mfilename,iface);
                                    end
                                end
                            case 2
                                %21
                                fprintf(fid,'%d %d %d     NLay NRow NCol 3D-placement\n',o(iGrp).placement);
                            otherwise
                                error('%s: illegal placementOption <<%d>> use 1..2',o(iGrp).placementOption);
                        end
                    end
                otherwise
                    error('%s: INPUTSTYLE %d is not supported <<%d>>',mfilename,INPUTSTYLE);
            end
            
            if ischar(fnameOrFid)
                fclose(fid);
            end
        end
        function o = getEndPoints(o,fname)
            % mpath_particleGroupObj.getEndPoints(endPontsFileName)
            % The endpoints generated by modpath will be loaded into the
            % respective mpath_particleGroupObj's
            %
            % TO 130219
            
            % First load endpoints totally
            endp = readEndPoints(fname);
            
            for iGrp = 1:numel(o)
                % Get the column header labels
                o(iGrp).endpColHdr = endp.colHdr;

                % Get the endpoints pertaining to each particle group object
                o(iGrp).endPoints = endp.P(endp.P(:,2)==iGrp,:);
            end
        end
        function o = getTsrPoints(o,fname)
            % mpath_particleGroupObj.getEndPoints(endPointsFileName)
            % The endpoints generated by modpath will be loaded into the
            % respective mpath_particleGroupObj's
            %
            % TO 130219
            
            % First load endpoints totally
            tsr = readTsrPoints(fname);
            
            for iGrp = 1:numel(o)
                % Get the column header labels
                o(iGrp).tsrColHdr = tsr.colHdr;

                % Get the endpoints pertaining to each particle group object
                o(iGrp).tsrPoints = tsr.P(tsr.P(:,5)==iGrp,:);
            end
        end
        function o = getPathLines(o,fname)
            % mpath_particleGroupObj.getPathLines(pathLinesFileName)
            % The pathLines generated by modpath will be loaded into the
            % respective mpath_particleGroupObj's
            %
            % TO 130219
            
            % First load endpoints totally
            
            fprintf('\n%s: Reading path lines ...\n',mfilename);
            
            pth = readPath(fname);
            
            fprintf('\nSorting pathline points ...');
            for iGrp = 1:numel(o)
                % Get the column header labels
                o(iGrp).pathColHdr = pth.colHdr;

                % Get the endpoints pertaining to each particle group object
                o(iGrp).pathPoints = pth.P(pth.P(:,2)==iGrp,:);
                o(iGrp).pathPoints = sortrows(o(iGrp).pathPoints,1);
            end
            fprintf(' done\n');
        end
        function endPointStatistics(o)
            % mpath_particleGroupObj.statistics;
            % provide statistics for particle groups.
            % Notice: first load endpoints into particle groups using
            % pgrp = pgrp.getEndPoints
            for iGrp = 1:numel(o)
                if isempty(o(iGrp).endPoints)
                    fprintf('Group nr %d name = %d has no endPoints.\n',iGrp,o(iGrp).name);
                    fprintf('apply method grp = grp.getEndpoints(endPointsFileName) first.\n');
                else
                    %3
                    fprintf('\n# Statistics for endpoints of group %d, ''%s''\n',iGrp,o(iGrp).name);
                    fprintf('Total nr of endpoints     = %4d\n',size(o(iGrp).endPoints,1));
                    fprintf('    pending               = %4d\n',sum(o(iGrp).endPoints(:,3)==0));
                    fprintf('    active                = %4d\n',sum(o(iGrp).endPoints(:,3)==1));
                    fprintf('    normally terminated   = %4d\n',sum(o(iGrp).endPoints(:,3)==2));
                    fprintf('    zoneTerminated        = %4d\n',sum(o(iGrp).endPoints(:,3)==3));
                    fprintf('    unreleased            = %4d\n',sum(o(iGrp).endPoints(:,3)==4));
                    fprintf('    stranded              = %4d\n',sum(o(iGrp).endPoints(:,3)==5));
                    fprintf('# ==============================\n');
                end
            end
        end
        function tsrStatistics(o)
            % mpath_particleGroupObj.tsrStatistics;
            % provide statistics for particle groups.
            % Notice: first load tsrPoints into particle groups using
            % pgrp = pgrp.getTsrPoints
            for iGrp = 1:numel(o)
                if isempty(o(iGrp).tsrPoints)
                    fprintf('Group nr %d name = %d has no tsrPointRecords.\n',iGrp,o(iGrp).name);
                    fprintf('apply method grp = grp.getTsrPoints(tsrPointsFileName) first.\n');
                else
                    %3
                    fprintf('\n# Statistics for tsr, group %d ''%s''\n',iGrp,o(iGrp).name);
                    fprintf('Total nr of points     = %d\n', size(o(iGrp).tsrPoints,1));
                    fprintf('Total nr of time steps = %d\n', max(o(iGrp).tsrPoints(:,2)));
                    fprintf('Min time               = %g\n', min(o(iGrp).tsrPoints(:,3)));
                    fprintf('Max time               = %g\n', max(o(iGrp).tsrPoints(:,3)));
                    fprintf('# ==============================\n');
                end
            end
        end

        function pathStatistics(o)
            %PATHSTATISTICS -- provide statistics for particle groups.
            %
            % USAGE:
            %    mpath_particleGroupObj.pathStatistics;
            %
            % provide statistics for particle groups.
            % Notice: first load pathPoints into particle groups using
            % pgrp = pgrp.getPathPoints
            for iGrp = 1:numel(o)
                if isempty(o(iGrp).pathPoints)
                    fprintf('Group nr %d name = %d has no pathLineRecords.\n',iGrp,o(iGrp).name);
                    fprintf('apply method grp = grp.getPathPoints(pathPointsFileName) first.\n');
                else
                    %3
                    fprintf('\n# Statistics for path file, group %d ''%s''\n',iGrp,o(iGrp).name);
                    fprintf('Total nr of points     = %d\n', size(o(iGrp).pathPoints,1));
                    fprintf('Total nr of time steps = %d\n', max(o(iGrp).pathPoints(:,4)));
                    fprintf('Min time               = %g\n', min(o(iGrp).pathPoints(:,5)));
                    fprintf('Max time               = %g\n', max(o(iGrp).pathPoints(:,5)));
                    fprintf('# ==============================\n');
                end
            end
        end
        function h = plotEndPoints(o,zone,lSpec)
            %PLOTENDPOINTS -- plots endpoints of ponts that started in given zones.
            %
            % USAGE:
            %    hdl = mpath_particleGroupObj.plotEndpoints([zone,[lSpec]])
            %
            % plots endpoints of ponts that started in given zones.
            % where lSpec is a lineSpec as valid for plot(x,y,lSpec)
            % plots the endpoints of the particle groups
            % Note first load endpoints with
            % pgrp = pgrp.getParticles
            %
            % TO 130219
            ix = strmatchi('xGF',o(1).endpColHdr);
            iy = strmatchi('yGF',o(1).endpColHdr);
            iz = strmatchi('zGF',o(1).endpColHdr);
            
            iZone = strmatchi('iZoneI',o(1).endpColHdr);
            if nargin<2, zone=[]; end
            if nargin<3
                if ischar(zone)
                    lSpec=zone;
                    zone =[];
                else
                    lSpec = '';
                end
            end
                
            
            for iGrp = numel(o):-1:1
                if isempty(zone)
                    I = true(size(o(iGrp).endPoints(:,1)));
                else
                    I = ismember(o(iGrp).endPoints(:,iZone),zone);
                end
                if ~isempty(I)
                    if isempty(lSpec)
                            lSpec = o(iGrp).lineSpec;
                    end
                    h(iGrp) = plot3(o(iGrp).endPoints(:,ix),...
                                    o(iGrp).endPoints(:,iy),...
                                    o(iGrp).endPoints(:,iz),lSpec);
                end
            end
        end
        function h = plotStartPoints(o,zone,lSpec)
            %PLOTSTARTPOINTS -- plot start poitns ending in specific zones
            %
            % USAGE:
            %    hdl = mpath_particleGroupObj.plotStartpoints([lSpec])
            %
            % plot starting locations of particles ending in given zone(s).
            % lSpec is a lineSpec as valid for plot(x,y,lSpec)
            % First load endpoints with
            % pgrp = pgrp.getParticles
            %
            % TO 130219
            
            ix = strmatchi('xGI',o(1).endpColHdr);
            iy = strmatchi('yGI',o(1).endpColHdr);
            iz = strmatchi('zGI',o(1).endpColHdr);
            
            iZone = strmatchi('iZoneF',o(1).endpColHdr);
            
            for iGrp = numel(o):-1:1
                if nargin<2 || ~exist('zone','var') || isempty(zone)
                    I = true(size(o(iGrp).endPoints(:,1)));
                else
                    I = ismember(o(iGrp).endPoints(:,iZone),zone);
                end
                if any(I)
                    if nargin<2, lSpec = o(iGrp).lineSpec; end
                    h(iGrp) = plot3(o(iGrp).endPoints(I,ix),...
                                    o(iGrp).endPoints(I,iy),...
                                    o(iGrp).endPoints(I,iz),lSpec);
                end
            end
        end
        
        
        function dispEndPoints(o)
            %DISPENDPOINTS -- prints the list of endpoints in a more readable fashion
            %
            % USAGE:
            %    mpath_particleGroupObj.dispEndPoints()
            %
            % prints the list of endpoints in a more readable fashion,
            % using formats pertaining to the individual columns.
            %
            % TO 130221
            
           % Column headers
           Hdr = {  'id','iGrp','iStatus',...
                    'tI','tF',...
                    'iGridI','iLayI','iRowI','iColI','iFaceI','iZoneI',...
                    'xLI','yLI','zLI','xGI','yGI','zGI',...
                    'iGridF','iLayF','iRowF','iColF','iFaceF','iZoneF',...
                    'xLF','yLF','zLF','xGF','yGF','zGF','label'};
           
           % Desirede format width
           Lfmt=[5 5 5 12 12,...
                 5 5 5 5 5 5,...
                 12 12 12 12 12 12,...
                 5 5 5 5 5 5,...
                 12 12 12 12 12 12];
            
            % Limit column header width to that of format-1
            for i=1:numel(Lfmt)
                Hdr{i} = Hdr{i}(1:min(Lfmt(i),numel(Hdr{i})));
            end
            
            % Format code (integer of floating point general)
            fmtcode = ['dddgg' 'dddddd' 'gggggg' 'dddddd' 'gggggg'];

            % Generate the format for printing and the overtall table title
            fmt='';
            ttl='';
            for i=1:numel(fmtcode)
                fmt = sprintf(sprintf('%%s %%%%%d%s',Lfmt(i),fmtcode(i)),fmt);
                ttl = sprintf(sprintf('%%s %%%ds',Lfmt(i)),ttl,Hdr{i});
            end
            
            % print title
            fprintf('%s\n',ttl);
                 
            % print table for all groups sequentially
            for iGrp = 1:numel(o)
%                 fprintf('Group "%s", nr %d, %d points\n',...
%                     o(iGrp).name,iGrp,size(o(iGrp).P,1));
                fprintf([fmt '\n'],o(iGrp).endPoints');
            end
        end

        
        
        function h = plotTsrPoints(o,timeNrs)
            % hdl = mpath_particleGroupObj.plotTsrPoints([lSpec])
            % where lSpec is a lineSpec as valid for plot(x,y,lSpec)
            % plots the time series points of the particle groups
            % Note first load time series points with
            % pgrp = pgrp.getParticles
            %
            % TO 130219
            
            ix = strmatchi('xG',o(1).tsrColHdr);
            iy = strmatchi('yG',o(1).tsrColHdr);
            iz = strmatchi('zG',o(1).tsrColHdr);
            
            if nargin>1
                timeNrs = round(sort(timeNrs(:)));
            end
            
            for iGrp = numel(o):-1:1
                mrk = mf_marker(iGrp);
                [~,m1] = unique(o(iGrp).tsrPoints(:,1),'first');
                [~,m2] = unique(o(iGrp).tsrPoints(:,1),'last');

                if nargin<2
                    timeNrs=1:numel(m1);
                else
                    timeNrs=timeNrs(timeNrs<numel(m1));
                end
                
                for it=timeNrs(:)'
                    if m2(it)==m1(it), continue; end
                    clr = mf_color(it);
                    if ~isempty(o(iGrp).tsrPoints)
                        h(iGrp) = plot3(o(iGrp).tsrPoints(m1(it):m2(it),ix),...
                                        o(iGrp).tsrPoints(m1(it):m2(it),iy),...
                                        o(iGrp).tsrPoints(m1(it):m2(it),iz),[mrk clr]);
                    end
                end
            end
            
        end
        function hdl = plotPath(o,lSpec)
            % hdl = mpath_particleGroupObj.plotPthPoints([lSpec])
            % where lSpec is a lineSpec as valid for plot(x,y,lSpec)
            % plots the time series points of the particle groups
            % Note first load time series points with
            % pgrp = pgrp.getParticles
            %
            % TO 130219
            
            fprintf('\nPlotting pathlines\n');
            ix = strmatchi('xG',o(1).pathColHdr);
            iy = strmatchi('yG',o(1).pathColHdr);
            iz = strmatchi('zG',o(1).pathColHdr);
            for iGrp = numel(o):-1:1
                fprintf('Group %2d, ',iGrp);
                if ~isempty(o(iGrp).pathPoints)
                    if nargin<2
                        lSpec = mf_color(iGrp,'brgkmc');
                    end
                    fprintf(' sorting,');
                    [~,m1] = unique(o(iGrp).pathPoints(:,1),'first');
                    [~,m2] = unique(o(iGrp).pathPoints(:,1),'last');
                    hdl(iGrp).h = NaN(size(m1));
                    fprintf(' plotting,');
                    for iL=1:length(m1)
                        hdl(iGrp).h(iL) = plot3(o(iGrp).pathPoints(m1(iL):m2(iL),ix),...
                                              o(iGrp).pathPoints(m1(iL):m2(iL),iy),...
                                              o(iGrp).pathPoints(m1(iL):m2(iL),iz),...
                                              lSpec);
                    end
                    fprintf(' done.\n');
                end
            end
        end
    end
end
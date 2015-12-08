classdef mpath_particleObj
    % Defines particles as an objects, one per particle.
    % You can have an aribitrary long vector of particle objects.
    properties
        grid = 1;
        groupNr
        releaseTime    % release time(s) for this point
        globalXYZ      % world xyz coordinates
        localXYZ       % local coordinates within cell
        LRC            % Layer Row and ColumnIndex
        trackingTime
        label='';

    end
    methods
        function o = mpath_particleObj(varargin)
            % generate particle objects for use in MPATH6
            % irrrespective of the inputStyle used.
            %
            % Same input as mapth_particleGroupObj, but here we explicitly
            % generate a list of individual particles.
            %
            % General input structure
            %  particleObj(gr ,zoneVals,varNm,value,varNm,value,...])
            %  zoneVals = {zoneDef,varNm,value,varNm,value,...; ...
            %              zoneDef,varNm,value,varNm,value,...; ...
            %              ...}
            %  zoneVals is a cell array in which each line defines a
            %  particleGroupObj.
            %  if zoneVals is a numeric array, then { } will be placed
            %  to join it with all remaining arguments into a cell array
            %  thus turning it into one obeying the previous definition.
            %  Hence
            %
            %  with one particleGroup, braces can be left out
            %  particleObj(gr,zoneDef,varNm,value,...)
            %  particleObj(gr,{zoneDef,varNm,value},...);
            %
            %  zoneDef can be
            %     1) a list of XYZ coordinates defining particle positions
            %        a mapth_particleGroupObj array
            %     2) a definition of a cell block [ix1 ix2 iy1 iy2 iz1 iz2]
            %        each definition may contain multiple lines to define a
            %        set of cell blocks.
            %     3) a logical zone array of size gr.size
            %     4) a logical zone array of size(gr.Ny,gr.Nx) but then the
            %        zoneDef{iGrp,2} must be the layer number
            %
            % USAGE1
            %  particles = mpath_particlesObj(gr,globalXYZ,varNm,value,varNm,value...});
            %  particles = mpath_particlesObj(gr,{globalXYZ,varNm,valye,varNm,value...});
            %
            %  globalXYZ is a 3-column array in which each line has [xW yW zW]
            %  All following inputs must be <<varNm,values>> pairs, where varNm
            %    is one of the fields of the mpath_particleObj (case insensitive,
            %    with autocompletion)
            %
            % USAGE2
            %  particles = mpath_particlesObj(gr,mpath_particlGroupObj,varNm,value,varNm,value...});
            %  particles = mpath_particlesObj(gr,{mpath_particleGroupObj,...,varNm,value,...},varNm,value...});
            %
            %  mpath_particleGroupObj is a set of particle objecs or particle group objects for
            %  which properties will be specified or respecified.
            %
            % USAGE4
            %  o = mpath_particlesObj(gr{{zoneArray3D,      ...,varNm,Value,...},varNm,value...)
            %  o = mpath_particlesObj(gr,{zoneArray2D,iLay, ...varNm,value,...},varNm,value...);
            %
            % USAGE5
            %  o = mpath_particlesObj(gr{,cellBlock,      ...,varNm,Value,...},varNm,value...)
            %  o = mpath_particlesObj(gr, cellBlock,      ...,varNm,Value,... ,varNm,value...)
            %   cellBlock = [ix1 ix2 iy1 iy2 iz1 iz2; ...]
            %
            %
            % EXAMPLE of a zoneVals cell array:
            %   ptcl = mpath_particleOjb(zoneArray,...
            %         'IFace',[2 4 5],'placement',[3 6],'groupNm',ducklake','releaseTime',[0 10 100]);
            %   ptcl = mpath_particleOjb(zoneArray,...
            %       {zoneArray1,zones1,'IFace',[2 4 5],'placement',[3 6],'groupNm',ducklake'   ,'releaseTime',[0 10 100];...
            %        zoneArray2,zones2,'IFace',[ 3 ]  ,'placement',[2 4],'groupNm','gooseField','releaseTime',[0 10 100]});
            %   ptcl = mpath_particleOjb([10 15 4 8 1 2],,...
            %         'IFace',[2 4 5],'placement',[3 6],'groupNm',ducklake','releaseTime',[0 10 100]);
            %
            %   Notice:
            %     The order of varNm,value pairs is immaterial.
            %    varNm names are not case sensitive and autocompletion ensures that the
            %    correct field is used when varNm uniquely defines this field.
            %       IFace,iface numbers to place particles within the cells of this group.
            %             0 = interior 3D, rest use 1..6.
            %       placement determines the distribution of particles across
            %           the cell faces (IFace>0) or within the cell (IFace==0)
            %                Always use a 3 element vector [nLay,nRow,nCol]
            %                or a scalar if nLay==nRow==nCol. Any 2 element
            %                vector will be expanded as follows
            %                [nLay,nRow]-->[nLay,nRow,nCol==nRow].
            %                nLay is the subdivision of the particle along
            %                the vertical coordinate of the cell. nRow is
            %                the subdivion of the particles along the row
            %                direction of the cell (y-direction) and nCol
            %                is the subdivision along the column direction
            %                of the cell (x-direction). The correct
            %                directions will be selected for each IFace.
            %                E.g. IFace==0 uses all three, IFace=1 or 2
            %                uses nLay,nRow, IFace=3 or 4 uses nLay,nCol and
            %                iFace=5 or 6 uses nRow,nCol.
            %                For instance 'placement',[7 7 7] is equivalent
            %                to 'placement',7 and 'placement',[7 7].
            %       Other parameters see fields of this mpath_particleObj.
            %
            % TO 130127
            
            if nargin<2, return; end
            
            [gr,varargin] = getType(varargin,'gridObj',[]);
            if isempty(gr)
                error('missing gridObj');
            end
            
            [zoneVals,varargin] = getType(varargin,'cell',{});
            if isempty(zoneVals)
                zoneVals = varargin(1);
                vararginToDo = varargin(2:end);
                %varargin=[];
            end

            % analyze input, generate mpath_particleObj
            % dummy objects, used as gtemplate to store global variable values
            % ptemplate is particleObj template
            % gtemplate is particleGroupObj template
            ptemplate = mpath_particleObj();
            ptemplate = setProps(ptemplate,vararginToDo{:});
            gtemplate = mpath_particleGroupObj();
            gtemplate = setProps(gtemplate,vararginToDo{:});

            NGrp = size(zoneVals,1);
            
            for iGrp = NGrp:-1:1
                zoneDef = zoneVals{iGrp,1};
                zVals   = zoneVals(iGrp,2:end);
                
                %[inputStyle,iLay,zoneRest,gridCellRegionOption] = getInputStyle(zoneDef,zoneVals(iGrp,2:end));

                switch class(zoneDef)
                
                case 'double'
                    if size(zoneDef,2)==3
                        % mpath_particlesObj(gr,globalXYZ[,varNm,value,varNm,value...]);
                        % if called with two arguments, i.e. numel(varargin)==1
                        % assume that varargin is a list of N 3D points in
                        % world coordinates: [xW(1:N) yW(1:N) zW(1:N)]

                        [ix,iy,iz,XL,YL,ZL,XW,YW,ZW] = xyzindex(zoneDef,gr);                    
                        for i=length(ix):-1:1
                            o(i)=ptemplate;
                            o(i).LRC      = [iz(i),iy(i),ix(i)];
                            o(i).localXYZ = [XL(i),YL(i),ZL(i)];
                            o(i).globalXYZ= [XW(i),YW(i),ZW(i)];
                        end
                        continue;
                    elseif size(zoneDef,2)==6
                        zoneDef(:,1:2) = min(sort(zoneDef(:,1:2),2),gr.Nx);
                        zoneDef(:,3:4) = min(sort(zoneDef(:,3:4),2),gr.Ny);
                        zoneDef(:,5:6) = min(sort(zoneDef(:,5:6),2),gr.Nlay);
                        zoneDef = max(zoneDef,1);
                        IDX = false(gr.size);
                        for j = 1:size(zoneDef,1)
                            IDX(zoneDef(j,3):zoneDef(j,4),...
                                zoneDef(j,1):zoneDef(j,2),...
                                zoneDef(j,5):zoneDef(j,6)) = true;
                        end
                        IDX = find(IDX);
                    else
                        error('Unknown option for zoneDef');
                    end
                            
                case 'logical'
                    % if called with three arguments, i.e. numel(varargin)==2
                    % assume that varargin{1}=zoneArray and varargin{2} is a cell array
                    % called zoneVals, which on each line (zoneVals(i,:)) defines the
                    % properties of a group of cells to be generated

                    if ~isempty(zVals)
                        switch class(zVals{1,1})
                            case 'char'
                                % gridCellRegionOption = 3;
                                % zoneDef is a 3D logial array of gr.size
                                IDX = find(zoneDef);
                            case 'double'
                                % gridCellRegionOption = 3;
                                % zoneDef must be a 2D mask of gr.size(1:2)
                                % zVals{1} then is the maskLayer number of
                                % this group.
                                zVals = zVals(2:end);
                                IDX = find(zoneDef(:,:,1)) + numel(zoneDef)*(iLay-1);
                            otherwise
                                error('unknown option zoneDef');
                        end
                    else
                        % zoneDef must be a 3D logical array of gr.size
                        IDX = find(zoneDef);
                    end
                    otherwise
                        error('unknown type option zVals');
                end
                    
                ptemplate           = setProps(ptemplate,zVals{:});

                % expand placement to guarantee it is a 3 element vector [nL,nR,nC]
                gtemplate.placement(end:3) = gtemplate.placement(end);

                %% Get the coordinates of the actual particles, using their placement

                % Count the number of particles to be place in each cell
                N=0;
                for iface = gtemplate.IFace
                    switch iface
                        case 0
                            npnts = prod(gtemplate.placement([1 2 3]));
                        case {1,2}
                            npnts = prod(gtemplate.placement([1 2]));
                        case {3,4}
                            npnts = prod(gtemplate.placement([1 3]));
                        case {5,6}
                            npnts = prod(gtemplate.placement([2 3]));
                        otherwise
                            error('%s: iface must be one of [0..6] not %d',iface);                                        
                    end
                    N = N + npnts;
                end
                    
                % numbe of cells with particles
                Ncells = numel(IDX);
                % Total number of particles is Ncells*N
                N = N * Ncells;

                o(N) = mpath_particleObj(); % preallocate
                IDX  = IDX(:)';
                N    = 0;
                % For the individual ifaces in the IFace vector
                for iface = gtemplate.IFace
                    [xR yR zR] = particleRelCoords(iface,gtemplate.placement);
                    Npnts = numel(xR);

                    for idx = IDX
                        xW  = gr.XM(idx)+gr.DX(idx)*xR;
                        yW  = gr.YM(idx)+gr.DY(idx)*yR;
                        zW  = gr.ZM(idx)+gr.DZ(idx)*zR;
                        lrc = cellIndices(idx,gr.size,'LRC');

                        for ip=1:Npnts
                            o(N+ip)          = ptemplate;
                            o(N+ip).globalXYZ= [xW(ip) yW(ip) zW(ip)];
                            o(N+ip).localXYZ = [xR(ip)+0.5 yR(ip)+0.5 zR(ip)+0.5];
                            o(N+ip).LRC      = lrc;
                        end
                        N=N+Npnts;
                    end                                
                end
            end
        end

        function writeLoc(o,fid)
            % mpath_particleObj.writeLoc(fid)
            %
            % write particle locations as required for inputStyle 2:
            % [gridNr,Lay,Row,Col,LocalX,LocalY,LocalZ,label]
            for i=1:numel(o)
                fprintf(fid,'%d %d %d %d %g %g %g %s\n',...
                    o(i).gridNr,o(i).LRC,o(i).localXYZ,o(i).label);
            end
        end
        function write(o,fid,ID,releaseTime,label)
            % mpath_particleObj.write(fid)
            %
            % write particles as required for inputStyle 1:
            % [ParticleID,groupNr,gridNr,Lay,Row,Col,LocalX,LocalY,LocalZ,releaseTime,label]
            for it=1:length(releastTime)
                for i=1:numel(o)
                    fprintf(fid,'%d %d %d %d %d %g %g %g %g %s\n',...
                        ID,o(i).gridNr,o(i).LRC,o(i).localXYZ,releaseTime(it),label);
                    ID=ID+1;
                end
            end
        end
    end
    
end
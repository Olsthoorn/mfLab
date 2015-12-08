classdef shapeObj
    %SHAPEOBJ -- read shapefiles and generates shapeObj
    %
    % USAGE: shp = shapeObj(shapefile,['-verbose')
    %   shapefile name does not require extensio
    %  '-v of -verbose' silences the function
    %
    % TO 140512
    
    properties
        Type
        TypeName
        Box
        NumParts
        NumPoints
        MRange
        MArray
        Parts
        PartTypes % MultiPatch
        Point     % PointZ,
        Points
        ZRange
        ZArray
        UserData
    end
    properties (Constant)
        %% Legal o types
        shapeTypes={
            0, 'Null Shape'
            1, 'Point'
            3, 'PolyLine'
            5, 'Polygon'
            8, 'MultiPoint'
            11,'PointZ'
            13,'PolyLineZ'
            15,'PolygonZ'
            18,'MultiPointZ'
            21,'PontM'
            23,'PolyLineM'
            25,'PolygonM'
            28,'MultiPointM'
            31,'MultiPatch'};

        %% PartTypes defined for MultiPatch
        partTypes={
            0,'TriangleStrip'
            1,'TriangleFan'
            2,'OuterRing'
            3,'InnerRing'
            4,'FirstRing'
            5,'Ring'};
    end
    methods
        function o = shapeObj(FName,varargin)
            %READSHP reads an ERSI o file
            %
            % Example:
            %    [o,SHAPE]=readshp(FName,['-verbose'])
            %
            %   option '-v' or '-verbose' silences the function
            %
            % ESRI o files are specified in
            % 'http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf'
            %
            % A shapefile stores nontopological geometry and attribtute  information
            % for the spatial features in a dataset. The geometry for a features is
            % stored as a o comprising a set of vector coordinates.
            % Because shapefiles do not have the processing overhead of a toplogical
            % data structure, the have advantages over other data sources such as
            % faster drawing speed and edit ability. Shapefiles handle single features
            % that overlap of that are noncontiguous.  They also typically require less
            % disk space and are easier to read or write (From ESRI shapefile.pdf).
            % Shapefiles can support point, line and area features. Area features are
            % stored as closed loop, double-digitized polygons. Attributes are held in
            % a dBase(R) format file. Each attribute record has a one-to-one
            % relationship with the associated o record.
            %
            % An ESRI shapefile consists of a main file (*.shp) an index file (*.shx)
            % and a dBASE file (*.dbf).
            % The main file is a direct-access, variable-record length
            % file in which each record describes a o with a list of its vertices.
            % In the index file, each record contains the offset of the corresponding
            % main file record from the beginning of the main file. The dBASE table
            % contains feature attributes with one record per feature. The one-to-one
            % relationship between geometry and attributes is based on record number.
            % Attribute records in the dBASE file must be in the same order as records
            % in the main file.
            % Example:
            %   basename.shp  % main file
            %   basename.shx  % index file
            %   basename.dbf  % dBASE table
            %
            % A shapefile stores integer (4 bytes) and double precision numbers (8
            % bytes). Positive and negative infinite and NaN are not allowed. Any
            % values smaller than 1e-38 is considered NaN by the shapefile reader.
            %
            % Main file consists of a file header followed by an arbitrary number of
            % [record header | record contents] sets.
            %
            % Shapefile contents can be diveded into
            %  1) Data related (Main file record contents, main file header's data
            %     description files (Shape, Type, Bounding Box etc).
            %     These fields are in LITTLE ENDIAN byte order.
            %  2) File management related (File and record lengths, record offsets, and so
            %     on. These fields are in BIG ENDIAN byte order.
            %
            %
            % o is a o object (or array of o objects)
            % SHAPE is overall info over all shapes and includes the dBASE info
            %
            % dbfread must be in the search path 
            % 
            % See also: writeSHP readExp writeExp plotShp
            %
            % TO 090729 130818 --> object

            % Copyright 2009 2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
            % under free software foundation GNU license version 3 or later
            
            if nargin==0
                return;
            end            
            
            headerLen = 100;
            
            [silence,varargin] = getWord(varargin,'-');
            verbose            = ~silence | getWord(varargin,'v');
            
            % reads o file Header and the dBase file
            SHAPE = readShapeFileHeader(FName,verbose);
            
            shxB  = SHAPE.shxB; % Pointer to open BIG ENDIAN index file
            shpL  = SHAPE.shpL; % Pointer to open LITTLE ENDIAN o file
            
            fseek(shxB,headerLen,'bof');  % Position pointer in shx index file
            
            o(SHAPE.NShapes) = shapeObj();
            
            for iSh=1:SHAPE.NShapes
                % reading offset of record in shp file from record in shx file, skip
                % rest of header in shx file
                shpOffset =fread(shxB,1,'int')*2;    fread(shxB,1,'int');

                % position pointer in shp file read using big endian read offset in shp
                % file and skip rest of record header
                fseek(shpL,shpOffset,'bof');    fread(shpL,2,'int');

                % reading record contents, using shpL Little endian pointer
                o(iSh).Type=fread(shpL,1,'int');
                j=find(o(iSh).Type==vertcat(o(iSh).shapeTypes{:,1}));
                if isempty(j), error('Illegal fileshapeType %d in file ''%s.shp''\n,',FName);
                else
                   o(iSh).TypeName=o(iSh).shapeTypes{j,2};
                end

               switch o(iSh).Type
                   case 0  % Point, skip, only Type was read
                   case 1  % Point
                      o(iSh).Points=fread(shpL,[1,2],'double');
                   case {3,5}  % 'PolyLine' 'Polygon'
                      o(iSh).Box       =fread(shpL,4,'double');
                      o(iSh).NumParts  =fread(shpL,1,'int');
                      o(iSh).NumPoints =fread(shpL,1,'int');
                      o(iSh).Parts     =fread(shpL,o(iSh).NumParts,'int');
                      o(iSh).Points    =fread(shpL,[2,o(iSh).NumPoints],'double')';
                   case 8, % 'MultiPoint'
                      o(iSh).Box       =fread(shpL,4,'double');
                      o(iSh).NumParts  =fread(shpL,1,'int');
                      o(iSh).NumPoints =fread(shpL,1,'int');
                      o(iSh).Points    =fread(shpL,[2,o(iSh).NumPoints],'double')';
                      o(iSh).MRange    =fread(shpL,[1,2],'double');
                      o(iSh).MArray    =fread(shpL,o(iSh).NumPoints,'double');
                   case 11, % 'PointZ'  [x y z Measure]
                      o(iSh).Point     =fread(shpL,[4,o(iSh).NumPoints],'double')';         
                   case {13,15} % 'PolyLineZ' 'PolygonZ'
                      o(iSh).Box       =fread(shpL,4,'double');
                      o(iSh).NumParts  =fread(shpL,1,'int');
                      o(iSh).NumPoints =fread(shpL,1,'int');
                      o(iSh).Parts     =fread(shpL,o(iSh).NumParts,'int')/2+1;
                      o(iSh).Points    =fread(shpL,[2,o(iSh).NumPoints],'double')';
                      o(iSh).ZRange    =fread(shpL,[1,2],'double');
                      o(iSh).ZArray    =fread(shpL,o(iSh).NumPoints,'double');
                      o(iSh).MRange    =fread(shpL,[1,2],'double');
                      o(iSh).MArray    =fread(shpL,o(iSh).NumPoints,'double');
                   case 18, % 'MultiPointZ'
                      o(iSh).Box       =fread(shpL,4,'double');
                      o(iSh).NumPoints =fread(shpL,1,'int');
                      o(iSh).Points    =fread(shpL,[2,o(iSh).NumPoints],'double')';
                      o(iSh).ZRange    =fread(shpL,[1,2],'double');
                      o(iSh).ZArray    =fread(shpL,o(iSh).NumPoints,'double');
                      o(iSh).MRange    =fread(shpL,[1,2],'double');
                      o(iSh).MArray    =fread(shpL,o(iSh).NumPoints,'double');
                   case 21, % 'PontM'          
                      o(iSh).Points    =fread(shpL,[1,2],'double')';
                      o(iSh).MArray    =fread(shpL,1,'double')';
                   case {23,25} % 'PolyLineM' and 'PolygonM'
                      o(iSh).Box       =fread(shpL,4,'double');
                      o(iSh).NumParts  =fread(shpL,1,'int');
                      o(iSh).NumPoints =fread(shpL,1,'int');
                      o(iSh).Parts     =fread(shpL,o(iSh).NumParts,'int')/2+1;
                      o(iSh).Points    =fread(shpL,[2,o(iSh).NumPoints],'double')';
                      o(iSh).MRange    =fread(shpL,[1,2],'double');
                      o(iSh).MArray    =fread(shpL,o(iSh).NumPoints,'double');
                   case 28, % 'MultiPointM'
                      o(iSh).Box       =fread(shpL,4,'double');
                      o(iSh).NumPoints =fread(shpL,1,'int');
                      o(iSh).Points    =fread(shpL,[2,o(iSh).NumPoints],'double')';
                      o(iSh).ZRange    =fread(shpL,[1,2],'double');
                      o(iSh).ZArray    =fread(shpL,o(iSh).NumPoints,'double');
                      o(iSh).MRange    =fread(shpL,[1,2],'double');
                      o(iSh).MArray    =fread(shpL,o(iSh).NumPoints,'double');
                   case 31, % 'MultiPatch'};
                      o(iSh).Box       =fread(shpL,4,'double');
                      o(iSh).NumParts  =fread(shpL,1,'int');
                      o(iSh).NumPoints =fread(shpL,1,'int');
                      o(iSh).Parts     =fread(shpL,o(iSh).NumParts,'int')/2+1;
                      o(iSh).PartTypes =fread(shpL,o(iSh).NumParts,'int');
                        for j=1:size(o.partTypes,1)
                            if ~any(o(iSh).PatTypes==o.partTypes(j,1))
                                error('Illegal partTypes found [outsed 0..5] in Multipatch, o Type 31\n');
                            end
                        end
                      o(iSh).Points    =fread(shpL,[2,o(iSh).NumPoints],'double')';
                      o(iSh).ZRange    =fread(shpL,[1,2],'double');
                      o(iSh).ZArray    =fread(shpL,o(iSh).NumPoints,'double');
                      o(iSh).MRange    =fread(shpL,[1,2],'double');
                      o(iSh).MArray    =fread(shpL,o(iSh).NumPoints,'double');
                   otherwise
                       fprintf('Non implemented o file %d skipped\n',o.Type)
               end

             %% DBF data (was read in readShapeFileHeader
                for iField=1:length(SHAPE.attributes)
                    if SHAPE.attributes(iField).fieldtype == 'N'
                        o(iSh).UserData.(SHAPE.attributes(iField).fieldname)=SHAPE.attributes(iField).values(iSh);
                    else
                        o(iSh).UserData.(SHAPE.attributes(iField).fieldname)=deblank(SHAPE.attributes(iField).values{iSh});
                    end
                end               
            end
            fclose(shpL);
            fclose(shxB);
        end
       
        function plotBox(o,lineSpec)
            %PLOTBOX: plot shapes see plotshp
            % example: shapeObj.plotBox(lineSpec);
            %
            if nargin<2
                lineSpec = 'r-';
            elseif ~isLineSpec(lineSpec)
                error('illegal lineSpec = %s\n',lineSpec);
            end
            for iSh=1:length(o)
                box=o(iSh).Box;
                plot(box([1 3 3 1 1]),box([2 2 4 4 2]),lineSpec); hold on
            end
        end
        
        function plot(o,varargin)
            %PLOTS: plot shapes see plotshp
            %
            % USAGE: shapeObj.plot(lineSpec,plotoptons);
            %
            % for plotoptions see attributes of plot
            
            if nargin<2
                lineSpec = '';
            else
                lineSpec = varargin{1};
                if ~isempty(varargin), varargin(1)=[]; end
            end
            if ~isLineSpec(lineSpec)
                error('illegal lineSpec = %s\n',lineSpec);
            end
            
            for iSh=1:length(o)
                switch o(1).Type
                    case 1
                      xy = vertcat(o.Points);
                      if isempty(lineSpec), lineSpec = 'ro'; end
                      plot(xy(:,1),xy(:,2),lineSpec,varargin{:});
                    otherwise
                        if isempty(lineSpec), lineSpec = 'b-'; end
                        for iParts=1:length(o(iSh).NumParts)
                            first=o(iSh).Parts(iParts)+1;
                            if iParts<o(iSh).NumParts
                                last=o(iSh).Parts(iParts+1);
                            else
                                last=o(iSh).NumPoints;
                            end
                            range=first:last;
                            plot(o(iSh).Points(range,1),o(iSh).Points(range,2),lineSpec,varargin{:});
                            hold on
                        end
                end
            end
        end
        

        function fill(o,varargin)
            %FILL: fill shapes see shape.plot shape.label
            %
            % USAGE: shapeObj.fill(plotoptons);
            %
            % for plotoptions see attributes of fill and patch
            % TO 140514
                        
            for iSh=1:length(o)
                switch o(1).Type
                    case 1
                        o.plot(varargin{:});
                        return;
                    otherwise
                        for iParts=1:length(o(iSh).NumParts)
                            first=o(iSh).Parts(iParts)+1;
                            if iParts<o(iSh).NumParts
                                last=o(iSh).Parts(iParts+1);
                            else
                                last=o(iSh).NumPoints;
                            end
                            range=first:last;
                            fill(o(iSh).Points(range,1),o(iSh).Points(range,2),varargin{:});
                            hold on
                        end
                end
            end
        end
        

        function label(o,FLDNM,varargin)
            %LABEL: plot shapes see plotshp
            %
            % USAGE: shapeObj.label(FLDNM[,xOff,yOff],plotOptions);
            %
            % for plotoptions see attributes of plot
            
            fields = fieldnames(o(1).UserData);
            iFld = strmatchi(FLDNM,fields,'exact');
            if  ~iFld(1), return; end
            
            FLDNM = fields{iFld};
            
            [xOff,varargin] = getNext(varargin,'double',0);
            [yOff,varargin] = getNext(varargin,'double',0);
            for io=1:length(o)
                value = o(io).UserData.(FLDNM);
                if isnumeric(value), fmt = '%g'; else fmt='%s'; end

                switch o(io).Type
                    case 1
                        text(o(io).Points(1)+xOff,o(io).Points(2)+xOff,...
                            sprintf(fmt,value),varargin{:});
                    otherwise
                        for iParts=1:length(o(io).NumParts)
                            first=o(io).Parts(iParts)+1;
                            if iParts<o(io).NumParts
                                last=o(io).Parts(iParts+1);
                            else
                                last=o(io).NumPoints;
                            end
                            range=first:last;
                            text(mean(o(io).Points(range,1))+xOff,...
                                 mean(o(io).Points(range,2))+yOff,...
                                 sprintf(fmt,value),varargin{:});
                            hold on
                        end
                end
            end
        end
        
        function table(o)
            %TABLE -- prints te attribute tabel for these shapes
            % USAGE: shape.table
            % TO 140513
            
            fprintf('\n');
            fields = fieldnames(o(1).UserData);
            for ifld=1:numel(fields)
                fprintf('\t%s',fields{ifld});
            end
            fprintf('\n');
            for io=1:numel(o)
                for ifld = 1:numel(fields)
                    if ischar(o(io).UserData.(fields{ifld}))
                        fprintf('\t%s',o(io).UserData.(fields{ifld}));
                    else
                        fprintf('\t%g',o(io).UserData.(fields{ifld}));
                    end
                end
                fprintf('\n');
            end
        end

 
       function o = setAttr(o,Key,hdrs,table,varargin)
            %SETATTR -- sets attributes for shapes
            %
            % USAGE: shape.setAttr(AttrNm,hdrs,table);
            %   replaces or addes the attribute values of the shapes
            %   the Key is the Attr name is must be in the UserData fields
            %   and in the hdrs (a cell array of attribute names). The
            %   table contains the values, the columns must match those of
            %   the headers.
            %
            %   New attribures can be added this way.
            %   Value of existing attributes can be changed this way
            %
            % TO 140513
            
            select = getWord(varargin,'select');
            
            fields = fieldnames(o(1).UserData);
            iU = strmatchi(Key,fields,'exact');
            iH = strmatchi(Key,hdrs ,'exact');
            
            
            if ~iU, error('%s not in User data of shapes',Key); end
            if ~iH, error('%s not in provided headers',   Key); end
            if numel(hdrs) ~= size(table,2)
                error('number of hdrs does not match number of columns of provided table');
            end
            
            % make sure other hdrs have correct spelling
            Key      = fields{iU};  % exactly same
            hdrs{iH} = Key; % exactly same spelling
            
            for ih=1:numel(hdrs)
                ifld = strmatchi(hdrs{ih},fields,'exact');
                if ifld, hdrs{ih} = fields{ifld}; end
            end
            
            if select
                for io=numel(o):-1:1
                    key(io) = o(io).UserData.(Key);
                end
                o = o(ismember(key,[table{:,iH}]));
            end
            
            % for all objects
            for io = 1:numel(o)
            % look up the line in the code table corresponding with the
               line = find(ismember(vertcat(table{:,iH}),o(io).UserData.(Key)));
               if isempty(line), continue; end
               for ih = 1:numel(hdrs)
                    if ih==iH, continue; end
                    o(io).UserData.(hdrs{ih}) = table{line,ih};
               end
            end
       end

      function stress = stress(o,type,gr)
            %STRESS -- sets attributes for shapes
            %
            % USAGE: STRESS = shape.stress(type,gr);
            %   replaces or addes the attribute values of the shapes
            %   the Key is the Attr name is must be in the UserData fields
            %   and in the hdrs (a cell array of attribute names). The
            %   table contains the values, the columns must match those of
            %   the headers.
            %
            %   New attribures can be added this way.
            %   Value of existing attributes can be changed this way
            %
            % TO 140513
            legalTypes  = {'GHB','DRN','RIV','CHD'};
            
            type = upper(type);
            if ~ismember(type,legalTypes)
                error('illegal type ''%s'', use oneof ''%s''',type,...
                    sprintf(' ''%s''',legalTypes{:}));
            end
            
            % in case Idx is missing, get it
            if ~o(1).UserData.Idx
                o = o.Idx(gr);
            end
            
            for io = numel(o):-1:1                
               u = ones(size(o(io).UserData.Idx));  
               I = o(io).UserData.Idx(:);
               
               % polyline have are lineOb use dL instead of Area
               if o(io).Type == 3 % polyline
                   dL = o(io).UserData.dL(:);
                   switch type
                        case 'GHB'
                               stress{io,1} = [I u*o(io).UserData.h dL./o(io).UserData.c ];               
                        case 'DRN'
                               stress{io,1} = [I u*o(io).UserData.h dL./o(io).UserData.c ];               
                        case 'RIV'
                               stress{io,1} = [I u*o(io).UserData.h dL./o(io).UserData.c u*o(io).UserData.hBot ];               
                        case 'CHD'
                               stress{io,1} = [I u*o(io).UserData.h([1 end])];               
                       otherwise
                   end                   
               else
                   % all others use Area
                    switch type
                        case 'GHB'
                               stress{io,1} = [I u*o(io).UserData.h gr.Area(I)/o(io).UserData.c ];               
                        case 'DRN'
                               stress{io,1} = [I u*o(io).UserData.h gr.Area(I)/o(io).UserData.c ];               
                        case 'RIV'
                               stress{io,1} = [I u*o(io).UserData.h gr.Area(I)/o(io).UserData.c u*o(io).UserData.hBot ];               
                        case 'CHD'
                               stress{io,1} = [I u*o(io).UserData.h([1 end])];               
                        otherwise
                            % skip
                    end
               end
            end
      end

      function  Qs = getQ(o,Q)
          %getQ -- retrieves injection into model from object
          % USAGE: Q = shp.setQ(Q)
          %   Q is from model, full cell array
          %   TO 140515
          
          o = o.setQ(Q);
          Qin = 0; 
          Qout= 0;
          Q   = 0; 
          for io=1:numel(o)
              Qin = Qin + o(io).UserData.Qin;
              Qout= Qout+ o(io).UserData.Qout;
              Q   = Q   + o(io).UserData.Q;
          end
          Qs= [Qin Qout Q];
      end

      function  [x,y] = points(o,n,gr)
          %points -- generates starting points for Moc
          % USAGE: [x,y] = shp.points(n)
          %   n = distributed points n x n per cell
          %  -n = n x n randomized points per cell
          %   TO 160515
          
        if nargin<3, error('not enough input arguments'); end
          
            randomize = n<0;
            n      = abs(n);
            Npnt   = n*n;
            u      = (1/(2*n):1/n:1-1/(2*n))-0.5;
            uu     = ones(n,1)  * u;
            vv     = u' * ones(1,n);
            for io=numel(o):-1:1
                Idx    = o(io).UserData.Idx;
                Ncells = numel(Idx);
                xm     = gr.Xm(Idx)*ones(1,Npnt);
                ym     = gr.Ym(Idx)*ones(1,Npnt);
                dx     = gr.dX(Idx)*ones(1,Npnt);
                dy     = gr.dY(Idx)*ones(1,Npnt);
                xr0    = ones(Ncells,1)*uu(:)';
                yr0    = ones(Ncells,1)*vv(:)';

                if randomize
                    xr0 = xr0 + rand(size(xr0))-0.5;
                    yr0 = yr0 + rand(size(yr0))-0.5;

                    %mirrors at cell faces to keep particles in their original cells
                    xr0(xr0<-0.5) = -1 - xr0(xr0<-0.5);
                    xr0(xr0> 0.5) =  1 - xr0(xr0> 0.5);
                    yr0(yr0<-0.5) = -1 - yr0(yr0<-0.5);
                    yr0(yr0> 0.5) =  1 - yr0(yr0> 0.5);
                end
                
                xy{io} = [xm(:) + dx(:).*xr0(:), ym(:) + dy(:).*yr0(:)];
            end
            xy = cell2list(xy);
            x  = xy(:,1); y=xy(:,2);
      end

      function  [x,y] = edgePoints(o,L)
          %edgePoints -- generates starting points for Moc on edges of shape
          % USAGE: [x,y] = shp.points(n)
          %   n = distributed points n x n per cell
          %  -n = n x n randomized points per cell
          %   TO 160515
          
        if nargin<2, error('not enough input arguments'); end
          
        for io=numel(o):-1:1
            x = o(io).Points(:,1);
            y = o(io).Points(:,2);

            s = cumsum([0; sqrt(diff(x).^2+diff(y).^2)]);
            sp= (s(1):L:s(end))';
            xy{io} = [interp1(s,x,sp) interp1(s,y,sp)];
        end
        xy = cell2list(xy);
        x  = xy(:,1); y=xy(:,2);
      end


      function N = captured(o,P)
          %CAPTURED -- count particles captured by shp
          % USAGE: N = shp.captured(P,gr)
          %    N(size t) partcles captured
          %    P = output of Moc (has P(it).x,P(it).y,P(it).Icells)
          N = zeros(1,numel(P));
          for it=numel(P):-1:1
              for io=1:numel(o)
                  N(it) = N(it) + sum(ismember(P(it).Icells,o(io).UserData.Idx));
              end
          end
      end
      
      function o = setQ(o,Q)
          %setQ -- inserts Q into object (infiltration positive)
          % USAGE shap = shape.setQ(Q)
          %  Q is from model, full cell array
          %  TO 140515
          
          for io=1:numel(o)
              dQ = Q(o(io).UserData.Idx);
              o(io).UserData.Q    = sum(dQ);
              o(io).UserData.Qin  = sum(dQ(dQ>0));
              o(io).UserData.Qout = sum(dQ(dQ<0));
          end
      end
       
       function [o,IDR] = Idx(o,gr)
           %Idx computes global grid indices whos cell centers fall inside
           %shape
           InAll = false(gr.size);
           for io=1:length(o)
               In    = false(gr.size);
                switch o(io).Type
                    case 1
                        o(io).UserData.Idx = xyzindex(o(io).Points(1),o(io).Points(2),gr.xGr,gr.yGr);
                    case 3 % line shape
                        % compute the distance to the line shape. Points
                        % closer than some threshold will be included
                        % The threshhold if cell size or absolute
                        % bufferen en punten opvragen voor elk lijnstuk
                        % cut line through grid
                        if ~isfield(o(io).UserData,'Idx')
                            o(io).UserData.Idx = [];
                            o(io).UserData.dL  = [];
                        end
                        
                        % compute cell face intersections using vector
                        % analysis
                        xL = o(io).Points(:,1);
                        yL = o(io).Points(:,2);
                        dx = diff(xL);
                        dy = diff(yL);
                        dL_ = sqrt(dx.^2+dy.^2);
                        for iL=numel(dL_):-1:1
                            % from current point to all grid lines
                            lambX = (gr.xGr - xL(iL)) ./ dx(iL);
                            lambY = (gr.yGr - yL(iL)) ./ dy(iL);
                            % ralative to vector length ([0 1]) adds start
                            % and end of vector itself.
                            lamb  = unique([0 1 lambX(:)' lambY(:)']);
                            % only within length of this vector
                            lamb  = lamb(lamb>=0 & lamb<=1);
                            % midpoints of line sections in cells
                            lambM = 0.5*(lamb(1:end-1) + lamb(2:end));
                            % midpoint coordinates
                            xm_   = xL(iL) + lambM * dx(iL);
                            ym_   = yL(iL) + lambM * dy(iL);
                            % global cell number of midpoints
                            IL{iL,1} = [xyzindex(xm_,ym_,gr.xGr,gr.yGr), ...
                                        diff(lamb)' * dL_(iL) ];
                        end
                        % make unique and add dL for repetited cell numbers
                        L          = sortrows(cell2list(IL));
                        [I,Ifirst] = unique(L(:,1),'first');
                        [~,Ilast ] = unique(L(:,1),'last');
                        o(io).UserData.Idx = I;
                        o(io).UserData.dL  = zeros(size(I));
                        for iRow=1:numel(I)
                            o(io).UserData.dL(iRow) = sum(L(Ifirst(iRow):Ilast(iRow),2));
                        end
                        
                    otherwise
                        for iParts=1:length(o(io).NumParts)
                            first=o(io).Parts(iParts)+1;
                            if iParts<o(io).NumParts
                                last=o(io).Parts(iParts+1);
                            else
                                last=o(io).NumPoints;
                            end
                            range=first:last;
                            In = In | inpolygon(gr.Xm,gr.Ym, ...
                                        o(io).Points(range,1),...
                                        o(io).Points(range,2));
                        end
                        o(io).UserData.Idx = find(In);
                end
                if o(io).Type ~= 1
                    InAll = InAll | In;
                elseif ~isnan(o(io).UserData.Idx)
                    InAll(o(io).UserData.Idx)=true;
                end
           end
           IDR = find(InAll);
       end

      function o = setHead(o,Phi)
           %SETHEAD -- Gets the ehead from Phi and stores it internally
           %
           % USAGE: shapeObj.setHead(Phi)
           %
           %  requires that shapObj knows its global cell number
           %   see shapeObj.Idx to get it
           %
           % it currently only works for point shapes
           % Phi currently is a 2D head array. Must be extended to 3D and
           % steady state. Workaround, call it repeatedly with the array
           % corresponding to a specific stress period.
           %
           % TO 140513

           for io=1:length(o)
                switch o(io).Type
                    case 1
                        o(io).UserData.Hd = Phi(o(io).UserData.Idx);
                    otherwise
                       % skip
                end
           end
      end

      function A = getA(o,A,FldNm)
           %SETHEAD -- gets the numeric array A to the value in atribute
           % FldNm  if FldNm is a string, or to FldNm if FldNm is numeric.
           % It does so for the intersectd cells, set with
           %  shapeObj = shapeObj.Idx(grdObj)
           %
           % USAGE: A      = shp.setA(A,FLDNM)
           %        IBOUND = shp.getA(IBOUND,-1);
           %        IH     = shp.getA(IH,'STAGE');
           %        C      = shp.getA(C, 'C');
           %
           %  requires that shapObj knows its global cell number
           %
           %  SEE ALSO shapeObj.Idx shapeObj.setHead
           %
           % it currently only works in 2D
           %
           % TO 140513
           
           if nargin<3,
               error('Insufficient input arguments, use IH = shp.getIh(IH,FldNm);');
           end
           if ~isnumeric(FldNm)
               fields = fieldnames(o(1).UserData);
               ih = strmatchi(FldNm,fields);
               if ~ih
                   error('shapeObj has no attribute ''%s'', use oneof ''%s''',...
                       FldNm,sprintf(' ''%s''',fields{:}));
               end          
               FldNm = fields{ih}; % exact spelling          
               for io=1:length(o)
                   if ~all(isnan(o(io).UserData.Idx) | isempty(o(io).UserData.Idx))
                       A(o(io).UserData.Idx) = o(io).UserData.(FldNm);
                   end
               end
           else
               Value = FldNm;
               for io=1:length(o)
                   if ~all(isnan(o(io).UserData.Idx) | isempty(o(io).UserData.Idx))
                       A(o(io).UserData.Idx) = Value;
                   end
               end
           end
      end


        function writeShp(o,FName)
            % WRITESHP Write ARCVIEW shapefile.
            % Writes data to ESRI shape file, generating three files
            %    FName.shx
            %    FName.shp
            %    FName.dbf
            % Standard fields in record:
            % [shapeType, x, y, ......] in lower case
            % Implemented types: point, polyline, en polygon.
            %
            % dBASE field names always in upper case
            %
            % Assumes dbfWrite in search path
            %
            % See also READSHP, READEXP, WRITEEXP, PLOTSHP
            % 
            % Frans Schaars, 2000

            ESRI.filecode=9994;
            ESRI.unused  =[0 0 0 0 0];
            ESRI.filelength=NaN; %weten we nog niet
            ESRI.version=1000;

            Zmin=0;
            Zmax=0;
            Mmin=0;
            Mmax=0;
            
            %% write big ENDIAN part of file header (see ESRI document)
            fidout=fopen([FName,'.shp'],'w','b'); % big ENDIAN

            fwrite(fidout,[ESRI.filecode ESRI.unused ESRI.filelength],'int');

            % write little ENDIAN part of file header (see ESRI document)
            fidout=fopen([FName,'.shp'],'a','l'); % little ENDIAN
            
            fseek( fidout,28,'bof');
            fwrite(fidout,ESRI.version,'int');

            switch o(1).shapeType
            case 'point'   
               fwrite(fidout,1,'int');
               Xmin=min([o.x]);
               Ymin=min([o.y]);
               Xmax=max([o.x]);
               Ymax=max([o.y]);

            case 'polyline'
               fwrite(fidout,3,'int');
               xtmp=(cat(2,o.x));
               ytmp=(cat(2,o.y));
               Xmin=min(cat(2,xtmp{:}));
               Ymin=min(cat(2,ytmp{:}));
               Xmax=max(cat(2,xtmp{:}));
               Ymax=max(cat(2,ytmp{:}));

            case 'polygon'
               fwrite(fidout,5,'int');
               xtmp=(cat(2,o.x));
               ytmp=(cat(2,o.y));
               Xmin=min(cat(2,xtmp{:}));
               Ymin=min(cat(2,ytmp{:}));
               Xmax=max(cat(2,xtmp{:}));
               Ymax=max(cat(2,ytmp{:}));

            end

            fwrite(fidout,[Xmin,Ymin,Xmax,Ymax,Zmin,Zmax,Mmin,Mmax],'double');

            o=50;
            b=108;
            for iSp=1:length(o)
               fwrite(fidout,[NaN NaN],'int');%hier komt later wat
               shapeType=o(iSp).shapeType;
               switch shapeType
               case 'point'
                  fwrite(fidout,1,'int');
                  fwrite(fidout,o(iSp).x,'double');
                  fwrite(fidout,o(iSp).y,'double');
                  co(iSp)=2+2*4;
               case 'polyline'
                  fwrite(fidout,3,'int');
                  fwrite(fidout,o(iSp).box,'double');
                  fwrite(fidout,o(iSp).numparts,'int');
                  fwrite(fidout,o(iSp).numpoints,'int');
                  fwrite(fidout,o(iSp).parts,'int');
                  XY=[];
                  for j=1:o(iSp).numparts
                     XY=[XY [o(iSp).x{j};o(iSp).y{j}]];
                  end
                  fwrite(fidout,XY,'double');
                  co(iSp)=22+o(iSp).numparts*2+o(iSp).numpoints*2*4;
               case 'polygon'
                  fwrite(fidout,5,'int');
                  fwrite(fidout,o(iSp).box,'double');
                  fwrite(fidout,o(iSp).numparts,'int');
                  fwrite(fidout,o(iSp).numpoints,'int');
                  fwrite(fidout,o(iSp).parts,'int');
                  XY=[];
                  for j=1:o(iSp).numparts
                     XY=[XY [o(iSp).x{j};o(iSp).y{j}]];
                  end
                  fwrite(fidout,XY,'double');
                  co(iSp)=22+o(iSp).numparts*2+o(iSp).numpoints*2*4;
               end
               beginrecord(iSp)=b;
               b=b+2*(co(iSp)+4);
               offset(iSp)=o;
               o=o+co(iSp)+4;
            end

            fidout=fopen([FName,'.shp'],'r+','b');
            for iSp=1:length(beginrecord)
               fseek(fidout,beginrecord(iSp)-8,'bof');
               fwrite(fidout,iSp,'int');
               fwrite(fidout,co(iSp),'int');
            end
            %filelength
            fseek(fidout,24,'bof');
            fwrite(fidout,sum(co+4)+50,'int');



            %SHX
            fidout=fopen([FName,'.shx'],'w','b');
            fwrite(fidout,filecode,'int');
            fwrite(fidout,unused,'int');
            fwrite(fidout,50+4*length(o),'int');

            fidout=fopen([FName,'.shx'],'a','l');
            fseek(fidout,28,'bof');
            fwrite(fidout,version,'int');

            switch o(1).shapeType
            case 'point'   
               fwrite(fidout,1,'int');
            case 'polyline'
               fwrite(fidout,3,'int');
            case 'polygon'
               fwrite(fidout,5,'int');
            end

            fwrite(fidout,[Xmin,Ymin,Xmax,Ymax,Zmin,Zmax,Mmin,Mmax],'double');

            fidout=fopen([FName,'.shx'],'a','b');

            for iSp=1:length(o)
               fwrite(fidout,offset(iSp),'int');
               fwrite(fidout,co(iSp),'int');
            end
            fclose all;

            %DBFfile
            fn=fieldnames(o);
            c=0;
            for iFld=1:length(fn)
               fldNm=fn{iFld};
               j = strmatchi(fldNm,fn);
               if upper(fldNm)==fldNm %hoofdletters
                  c=c+1;
                  eval(['tmp=o(1).',fldNm,';'])
                  if isnumeric(tmp)
                     eval(['mat=cat(1,o.',fldNm,');'])
                     Type='N';
                     cw=16;
                     dec=2;
                  else
                     eval(['mat=char(o.',fldNm,');'])
                     Type='C';
                     cw=size(mat,2);
                     dec=0;
                  end
                  data{1,c}=fldNm;
                  data{2,c}=Type;
                  data{3,c}=cw;
                  data{4,c}=dec;
                  data{5,c}=mat;
               end
            end   
            dbfwrite([FName,'.dbf'],data);
        end
    end
end

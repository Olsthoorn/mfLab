classdef esriShapeObj
    %ESRISHAPEOBJ class def for ESRI shape objects
    %
    properties
        type
        typeName = 'unknown';
        parts
        points
        zArray
        mArray
        created
        filename = 'unknown';
    end
    properties (Dependent = true)
    end
    methods
        function o = esriShapeObj(shapefileNm)
            %ESRISHAPEOBJ constructor to generate an ESRI shape object
            %
            o.created = now;
            switch nargin
                case 0, return;
                case 1,
                    % only if nargin==1, we expect the argument to be the
                    % name of a shape file
                    o.filename = shapefileNm;
                    if exist(o.filename,'file')
                        o = o.readShp(o.filename);
                    else
                        error('%s: Can''t open file <<%s>> to read objects of class <<%s>>',...
                            mfilename,o.filename,class(o));
                    end
                otherwise
                    return;  % read from worksheet rather than shape file
            end
        end
        
        function [o,data]=readShp(o,filename)
        % [o,SHAPE]=readshp(filename) ---- use no extension
        % Reads ESRI o file, specified in
        % 'http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf'
        % A shapefile consists of a mainfile an strmatchi file and a dbase file
        % e.g. basename.shp, basename.shx, basename.dbf, basename the same, fnames 8.3 convention
        %
        % o is a struct array of shapes which includes the o data extracted
        % from the dbf file
        % SHAPE is overall info over all shapes
        % data is the contents of the dbf file (not needed, also per record in o-struct)
        %
        % dbfread must be in the search path 
        % 
        % See also WRITESHP, READEXP, WRITEEXP, PLOTSHP
        %
        % TO 090729


        % Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
        % under free software foundation GNU license version 3 or later

        %%
        fprintf('\n\nreadshp --- %s\n',datestr(now));

        [~,filename,~] =fileparts(filename); %  cutoff file extension

        o.filename = filename;
        o.created = now;

        headerlen=100;  % for shp and shx files

        %% Legal o types
        shapetypes={
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
        parttypes={
            0,'TriangleStrip'
            1,'TriangleFan'
            2,'OuterRing'
            3,'InnerRing'
            4,'FirstRing'
            5,'Ring'};

        %% Skip the header because it is idential to that of the .shp file, where
        %% we'll read it

        shxB =fopen([filename,'.shx'],'r','b');  % big endian bite ordering (Sun or Motorola)
        shpB =fopen([filename,'.shp'],'r','b');  % big endian bite ordering (Sun or Motorola)
        shpL =fopen([filename,'.shp'],'r','l');  % little endian (Intel or PC)

        if shxB<1, error('Can''t open file ''%s.shx''\n',filename); end
        if shpB<1, error('Can''t open file ''%s.shp''\n',filename); end
        if shpL<1, error('Can''t open file ''%s.shp''\n',filename); end

        shxHdr=fread(shxB,100,'int8');
        shpHdr=fread(shpB,100,'int8');
        if ~all(shxHdr([1:25 29:100])==shpHdr([1:25 29:100])), error('Header in .shx and .shp file differ, must be the same\n'); end

        %% Header shx file
        fseek(shxB, 0,'bof'); fprintf('FileCode  = %d\n',fread(shxB,1,'uint'));
        fseek(shxB,24,'bof'); fprintf('FileLenShx= %d\n',fread(shxB,1, 'int'));  % number of bytes
        fseek(shxB,28,'bof'); fprintf('Version   = %d\n',fread(shxB,1, 'int'));

        %% Header shp file
        fseek(shpB, 0,'bof'); fprintf('FileCode  = %d\n',fread(shpB,1,'uint'));
        fseek(shpB,24,'bof'); fprintf('FileLenSp = %d\n',fread(shpB,1, 'int')); % in 16 bit words
        fseek(shpB,28,'bof'); fprintf('Version   = %d\n',fread(shpB,1, 'int'));

        fclose(shpB); % no longer needed

        %% Get number of o records from strmatchi file
        fseek(shxB,0,'eof');
        NShapes=(ftell(shxB)-headerlen)/8;

        %% read overall bounding box
%         fseek(shpL,36,'bof');
%             BB.Xmin=fread(shpL,1,'double');
%             BB.Ymin=fread(shpL,1,'double');
%             BB.Xmax=fread(shpL,1,'double');
%             BB.Ymax=fread(shpL,1,'double');
%             BB.Zmin=fread(shpL,1,'double');
%             BB.Zmax=fread(shpL,1,'double');
%             BB.Mmin=fread(shpL,1,'double');
%             BB.Mmax=fread(shpL,1,'double');

        %% Reading data

        fseek(shxB,headerlen,'bof');  % Position pointer in shx strmatchi file
        
        for iShape=NShapes:-1:1
            % reading offset of record in shp file from record in shx file, skip
            % rest of header in skx file
            shpOffset =fread(shxB,1,'int')*2;    fread(shxB,1,'int');

            % position pointer in shp file read using big endian read offset in shp
            % file and skip rest of record header
            fseek(shpL,shpOffset,'bof');    fread(shpL,2,'int');

            % reading record contents, using shpL Little endian pointer
            o(iShape).type=fread(shpL,1,'int');            
            
            j=find(o(iShape).type==vertcat(shapetypes{:,1}));
            
            if isempty(j), error('Illegal fileshapetype %d in file ''%s.shp''\n,',filename);
            else
               o(iShape).typeName=shapetypes{j,2};
            end

           switch o(iShape).type
               case 0  % Point, skip, only type was read
               case 1  % Point
                  o(iShape).points=fread(shpL,[1,2],'double');
               case {3,5}  % 'PolyLine' 'Polygon'
                  o(iShape).box      =fread(shpL,4,'double');
                  o(iShape).numParts =fread(shpL,1,'int');
                  o(iShape).numPoints=fread(shpL,1,'int');
                  o(iShape).parts    =fread(shpL,o(iShape).numParts,'int');
                  o(iShape).points   =fread(shpL,[2,o(iShape).numPoints],'double')';
               case 8, % 'MultiPoint'
                  o(iShape).box         =fread(shpL,4,'double');
                  o(iShape).numParts    =fread(shpL,1,'int');
                  o(iShape).numPoints   =fread(shpL,1,'int');
                  o(iShape).points      =fread(shpL,[2,o(iShape).numPoints],'double')';
                  o(iShape).mRange      =fread(shpL,[1,2],'double');
                  o(iShape).mArray      =fread(shpL,o(iShape).numPoints,'double');
               case 11, % 'PointZ'  [x y z Measure]
                  o(iShape).Point     =fread(shpL,[4,o(iShape).numPoints],'double')';         
               case {13,15} % 'PolyLineZ' 'PolygonZ'
                  o(iShape).box      =fread(shpL,4,'double');
                  o(iShape).numParts =fread(shpL,1,'int');
                  o(iShape).numPoints=fread(shpL,1,'int');
                  o(iShape).parts    =fread(shpL,o(iShape).numParts,'int')/2+1;
                  o(iShape).points   =fread(shpL,[2,o(iShape).numPoints],'double')';
                  o(iShape).zRange    =fread(shpL,[1,2],'double');
                  o(iShape).zArray    =fread(shpL,o(iShape).numPoints,'double');
                  o(iShape).mRange    =fread(shpL,[1,2],'double');
                  o(iShape).mArray    =fread(shpL,o(iShape).numPoints,'double');
               case 18, % 'MultiPointZ'
                  o(iShape).box      =fread(shpL,4,'double');
                  o(iShape).numPoints=fread(shpL,1,'int');
                  o(iShape).points   =fread(shpL,[2,o(iShape).numPoints],'double')';
                  o(iShape).zRange    =fread(shpL,[1,2],'double');
                  o(iShape).zArray    =fread(shpL,o(iShape).numPoints,'double');
                  o(iShape).mRange    =fread(shpL,[1,2],'double');
                  o(iShape).mArray    =fread(shpL,o(iShape).numPoints,'double');
               case 21, % 'PontM'          
                  o(iShape).points   =fread(shpL,[1,2],'double')';
                  o(iShape).mArray   =fread(shpL,1,'double')';
               case {23,25} % 'PolyLineM' and 'PolygonM'
                  o(iShape).box      =fread(shpL,4,'double');
                  o(iShape).numParts =fread(shpL,1,'int');
                  o(iShape).numPoints=fread(shpL,1,'int');
                  o(iShape).parts    =fread(shpL,o(iShape).numParts,'int')/2+1;
                  o(iShape).points   =fread(shpL,[2,o(iShape).numPoints],'double')';
                  o(iShape).mRange      =fread(shpL,[1,2],'double');
                  o(iShape).mArray      =fread(shpL,o(iShape).numPoints,'double');
               case 28, % 'MultiPointM'
                  o(iShape).box       =fread(shpL,4,'double');
                  o(iShape).numPoints =fread(shpL,1,'int');
                  o(iShape).points    =fread(shpL,[2,o(iShape).numPoints],'double')';
                  o(iShape).zRange    =fread(shpL,[1,2],'double');
                  o(iShape).zArray    =fread(shpL,o(iShape).numPoints,'double');
                  o(iShape).mRange    =fread(shpL,[1,2],'double');
                  o(iShape).mArray    =fread(shpL,o(iShape).numPoints,'double');
               case 31, % 'MultiPatch'};
                  o(iShape).box      =fread(shpL,4,'double');
                  o(iShape).numParts =fread(shpL,1,'int');
                  o(iShape).numPoints=fread(shpL,1,'int');
                  o(iShape).parts    =fread(shpL,o(iShape).numParts,'int')/2+1;
                  o(iShape).PartTypes=fread(shpL,o(iShape).numParts,'int');
                    for j=1:size(parttypes,1)
                        if ~any(o(iShape).PatTypes==parttypes(j,1))
                            error('Illegal parttypes found [outsed 0..5] in Multipatch, o type 31\n');
                        end
                    end
                  o(iShape).points   =fread(shpL,[2,o(iShape).numPoints],'double')';
                  o(iShape).zRange    =fread(shpL,[1,2],'double');
                  o(iShape).zArray    =fread(shpL,o(iShape).numPoints,'double');
                  o(iShape).mRange    =fread(shpL,[1,2],'double');
                  o(iShape).mArray    =fread(shpL,o(iShape).numPoints,'double');
               otherwise
                   fprintf('Non implemented o file %d skipped\n',o.type)
           end
        end
        fclose(shpL);

        %plot shapes see plotshp
        % for iShape=1:length(o)
        %     for iParts=1:length(o(iShape).numParts)
        %         first=o(iShape).parts(iParts)+1;
        %         if iParts<o(iShape).numParts
        %             last=o(iShape).parts(iParts+1);
        %         else
        %             last=o(iShape).numPoints;
        %         end
        %         range=first:last;
        %         box=o(iShape).box;
        %         plot(box([1 3 3 1 1]),box([2 2 4 4 2]),'r'); hold on
        %         plot(o(iShape).points(range,1),o(iShape).points(range,2));
        %         hold on
        %     end
        % end

        %% DBF data
        data=dbfread(filename)'; %PK Transformatie extra toegepast

        for iShape=NShapes:-1:1
            for iField=length(data):-1:1
                        o(iShape).(data(iField).fieldname)=data(iField).values(iShape);
            end
        end
        end
    end
end
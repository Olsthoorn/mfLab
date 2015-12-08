function [shape,SHAPE,data]=readshp(varargin)
%READSHP reads an ERSI shape file
%
% Example:
%    [shape,SHAPE,DBF]=readshp(FName)
%    [shape,SHAPE,DBF]=readshp(FName,'-v'); % verbose
%
% ESRI shape files are specified in
% 'http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf'
%
% A shapefile stores nontopological geometry and attribtute  information
% for the spatial features in a dataset. The geometry for a features is
% stored as a shape comprising a set of vector coordinates.
% Because shapefiles do not have the processing overhead of a toplogical
% data structure, the have advantages over other data sources such as
% faster drawing speed and edit ability. Shapefiles handle single features
% that overlap of that are noncontiguous.  They also typically require less
% disk space and are easier to read or write (From ESRI shapefile.pdf).
% Shapefiles can support point, line and area features. Area features are
% stored as closed loop, double-digitized polygons. Attributes are held in
% a dBase(R) format file. Each attribute record has a one-to-one
% relationship with the associated shape record.
%
% An ESRI shapefile consists of a main file (*.shp) an index file (*.shx)
% and a dBASE file (*.dbf).
% The main file is a direct-access, variable-record length
% file in which each record describes a shape with a list of its vertices.
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
% shape is a struct array of shapes which includes the shape data extracted
% from the dbf file
% SHAPE is overall info over all shapes
% data is the contents of the dbf file (not needed, also per record in shape-struct)
%
% dbfread must be in the search path 
% 
% See also: writeSHP readExp writeExp plotShp
%
% TO 090729 151129
%
% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%%

[FName, varargin] = getType(varargin,'char',[]);
[verbose, ~     ] = getType(varargin,'char',false);

if verbose, fprintf('\n\nreadshp --- %s\n',datestr(now)); end

if isempty(FName), error('name of shape file not given'); end

[P,FName ] =fileparts(FName);
FName = fullfile(P,FName);

headerlen=100;  % for shp and shx files

%% Legal shape types
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


%% Open file pointers using big and little endian as needed

% big endian file pointer into index file
shxB =fopen([FName,'.shx'],'r','b');  % big endian bite ordering (Sun or Motorola)

% big endian file pointer into main file
shpB =fopen([FName,'.shp'],'r','b');  % big endian bite ordering (Sun or Motorola)

% little endian time pointer into main file
shpL =fopen([FName,'.shp'],'r','l');  % little endian (Intel or PC)

if shxB<1, error('Can''t open file ''%s.shx''\n',FName); end
if shpB<1, error('Can''t open file ''%s.shp''\n',FName); end
if shpL<1, error('Can''t open file ''%s.shp''\n',FName); end

%% Skip the header because it is idential to that of the .shp file

shxHdr=fread(shxB,100,'int8');
shpHdr=fread(shpB,100,'int8');

% assert that the two headers are equal
if ~all(shxHdr([1:25 29:100])==shpHdr([1:25 29:100])),
    error('Header in .shx and .shp file differ, must be the same\n');
end

%% Header shx file
fseek(shxB, 0,'bof'); if verbose, fprintf('FileCode  = %d\n',fread(shxB,1,'uint')); end
fseek(shxB,24,'bof'); if verbose, fprintf('FileLenShx= %d\n',fread(shxB,1, 'int')); end % number of bytes
fseek(shxB,28,'bof'); if verbose, fprintf('Version   = %d\n',fread(shxB,1, 'int')); end

%% Header shp file
fseek(shpB, 0,'bof'); if verbose, fprintf('FileCode  = %d\n',fread(shpB,1,'uint')); end
fseek(shpB,24,'bof'); if verbose, fprintf('FileLenShp= %d\n',fread(shpB,1, 'int')); end % in 16 bit words
fseek(shpB,28,'bof'); if verbose, fprintf('Version   = %d\n',fread(shpB,1, 'int')); end

fclose(shpB); % big endian pointer into shape file no longer needed

%% Get number of shape records from index file

SHAPE.name=FName;

fseek(shxB,0,'eof');
SHAPE.NShapes=(ftell(shxB)-headerlen)/8;

%% read bounding box from shapefile header
fseek(shpL,32,'bof');
    SHAPE.type=fread(shpL,1,'int');
    if verbose, fprintf('Shape type = %s\n',SHAPE.type); end
    
fseek(shpL,36,'bof');
    SHAPE.Xmin=fread(shpL,1,'double');
    SHAPE.Ymin=fread(shpL,1,'double');
    SHAPE.Xmax=fread(shpL,1,'double');
    SHAPE.Ymax=fread(shpL,1,'double');
    SHAPE.Zmin=fread(shpL,1,'double');
    SHAPE.Zmax=fread(shpL,1,'double');
    SHAPE.Mmin=fread(shpL,1,'double');
    SHAPE.Mmax=fread(shpL,1,'double');

%% Reading data

fseek(shxB,headerlen,'bof');  % Position pointer in shx index file
shape(SHAPE.NShapes).Type=0; % allocate
for iShape=1:SHAPE.NShapes
    % reading offset of record in shp file from record in shx file, skip
    % rest of header in skx file
    shpOffset =fread(shxB,1,'int')*2;    fread(shxB,1,'int');
   
    % position pointer in shp file read using big endian read offset in shp
    % file and skip rest of record header
    fseek(shpL,shpOffset,'bof');    fread(shpL,2,'int');
    
    % reading record contents, using shpL Little endian pointer
    shape(iShape).Type=fread(shpL,1,'int');
    j=find(shape(iShape).Type==vertcat(shapetypes{:,1}));
    if isempty(j), error('Illegal fileshapetype %d in file ''%s.shp''\n,',FName);
    else
       shape(iShape).TypeName=shapetypes{j,2};
    end
  
   switch shape(iShape).Type
       case 0  % Point, skip, only type was read
       case 1  % Point
          shape(iShape).Points=fread(shpL,[1,2],'double');
       case {3,5}  % 'PolyLine' 'Polygon'
          shape(iShape).Box       =fread(shpL,4,'double');
          shape(iShape).NumParts  =fread(shpL,1,'int');
          shape(iShape).NumPoints =fread(shpL,1,'int');
          shape(iShape).Parts     =fread(shpL,shape(iShape).NumParts,'int');
          shape(iShape).Points    =fread(shpL,[2,shape(iShape).NumPoints],'double')';
       case 8, % 'MultiPoint'
          shape(iShape).Box       =fread(shpL,4,'double');
          shape(iShape).NumParts  =fread(shpL,1,'int');
          shape(iShape).NumPoints =fread(shpL,1,'int');
          shape(iShape).Points    =fread(shpL,[2,shape(iShape).NumPoints],'double')';
          shape(iShape).MRange    =fread(shpL,[1,2],'double');
          shape(iShape).MArray    =fread(shpL,shape(iShape).NumPoints,'double');
       case 11, % 'PointZ'  [x y z Measure]
          shape(iShape).Point     =fread(shpL,[4,shape(iShape).NumPoints],'double')';         
       case {13,15} % 'PolyLineZ' 'PolygonZ'
          shape(iShape).Box       =fread(shpL,4,'double');
          shape(iShape).NumParts  =fread(shpL,1,'int');
          shape(iShape).NumPoints =fread(shpL,1,'int');
          shape(iShape).Parts     =fread(shpL,shape(iShape).NumParts,'int')/2+1;
          shape(iShape).Points    =fread(shpL,[2,shape(iShape).NumPoints],'double')';
          shape(iShape).ZRange    =fread(shpL,[1,2],'double');
          shape(iShape).ZArray    =fread(shpL,shape(iShape).NumPoints,'double');
          shape(iShape).MRange    =fread(shpL,[1,2],'double');
          shape(iShape).MArray    =fread(shpL,shape(iShape).NumPoints,'double');
       case 18, % 'MultiPointZ'
          shape(iShape).Box       =fread(shpL,4,'double');
          shape(iShape).NumPoints =fread(shpL,1,'int');
          shape(iShape).Points    =fread(shpL,[2,shape(iShape).NumPoints],'double')';
          shape(iShape).ZRange    =fread(shpL,[1,2],'double');
          shape(iShape).ZArray    =fread(shpL,shape(iShape).NumPoints,'double');
          shape(iShape).MRange    =fread(shpL,[1,2],'double');
          shape(iShape).MArray    =fread(shpL,shape(iShape).NumPoints,'double');
       case 21, % 'PontM'          
          shape(iShape).Points    =fread(shpL,[1,2],'double')';
          shape(iShape).MArray    =fread(shpL,1,'double')';
       case {23,25} % 'PolyLineM' and 'PolygonM'
          shape(iShape).Box       =fread(shpL,4,'double');
          shape(iShape).NumParts  =fread(shpL,1,'int');
          shape(iShape).NumPoints =fread(shpL,1,'int');
          shape(iShape).Parts     =fread(shpL,shape(iShape).NumParts,'int')/2+1;
          shape(iShape).Points    =fread(shpL,[2,shape(iShape).NumPoints],'double')';
          shape(iShape).MRange    =fread(shpL,[1,2],'double');
          shape(iShape).MArray    =fread(shpL,shape(iShape).NumPoints,'double');
       case 28, % 'MultiPointM'
          shape(iShape).Box       =fread(shpL,4,'double');
          shape(iShape).NumPoints =fread(shpL,1,'int');
          shape(iShape).Points    =fread(shpL,[2,shape(iShape).NumPoints],'double')';
          shape(iShape).ZRange    =fread(shpL,[1,2],'double');
          shape(iShape).ZArray    =fread(shpL,shape(iShape).NumPoints,'double');
          shape(iShape).MRange    =fread(shpL,[1,2],'double');
          shape(iShape).MArray    =fread(shpL,shape(iShape).NumPoints,'double');
       case 31, % 'MultiPatch'};
          shape(iShape).Box       =fread(shpL,4,'double');
          shape(iShape).NumParts  =fread(shpL,1,'int');
          shape(iShape).NumPoints =fread(shpL,1,'int');
          shape(iShape).Parts     =fread(shpL,shape(iShape).NumParts,'int')/2+1;
          shape(iShape).PartTypes =fread(shpL,shape(iShape).NumParts,'int');
            for j=1:size(parttypes,1)
                if ~any(shape(iShape).PatTypes==parttypes(j,1))
                    error('Illegal parttypes found [outsed 0..5] in Multipatch, shape type 31\n');
                end
            end
          shape(iShape).Points    =fread(shpL,[2,shape(iShape).NumPoints],'double')';
          shape(iShape).ZRange    =fread(shpL,[1,2],'double');
          shape(iShape).ZArray    =fread(shpL,shape(iShape).NumPoints,'double');
          shape(iShape).MRange    =fread(shpL,[1,2],'double');
          shape(iShape).MArray    =fread(shpL,shape(iShape).NumPoints,'double');
       otherwise
           fprintf('Non implemented shape file %d skipped\n',shape.type)
   end
end
fclose(shpL);

%% DBF data
data=dbfread(FName)'; %PK Transformatie extra toegepast

for iShape=1:SHAPE.NShapes
    for iField=1:length(data)
                dataItem = data(iField).values(iShape);
                if isa(dataItem,'cell')
                    dataItem = dataItem{1};
                    if isa(dataItem,'char')
                        dataItem = strtrim(dataItem);
                    end
                end
                shape(iShape).Attributes.(data(iField).fieldname)=dataItem;
    end
end

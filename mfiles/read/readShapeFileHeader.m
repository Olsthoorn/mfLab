function SHAPE = readShapeFileHeader(FName,verbose)
%%READSHAPEFILEHEADER -- reads the header of an ESRI shapefile
%
% USAGE: readShapeFileHeader(FName,['-verbose'])
%
% called by shapeObj
%
% TO 140512

    if verbose, fprintf('\n\nreadshp --- %s\n',datestr(now)); end

    if strcmpi(FName(end-3:end),'.shp'), FName=FName(1:end-4); end  % cut off .shp if necessary

    headerLen=100;  % for shp and shx files

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

    shxHdr=fread(shxB,headerLen,'int8');
    shpHdr=fread(shpB,headerLen,'int8');

    % assert that the two headers are equal
    if ~all(shxHdr([1:25 29:headerLen])==shpHdr([1:25 29:headerLen])),
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

    fseek(shxB,0,'eof'); % got eof
    SHAPE.NShapes = (ftell(shxB)-headerLen)/8;

    fclose(shpB); % big endian pointer into o file no longer needed

    %% Get number of o records from index file

    SHAPE.Name=FName;
    
    %% read bounding box from shapefile header
    fseek(shpL,32,'bof');
        SHAPE.Type=fread(shpL,1,'int');
        if verbose, fprintf('Shape type = %s\n',SHAPE.Type); end

    fseek(shpL,36,'bof');
        SHAPE.Xmin=fread(shpL,1,'double');
        SHAPE.Ymin=fread(shpL,1,'double');
        SHAPE.Xmax=fread(shpL,1,'double');
        SHAPE.Ymax=fread(shpL,1,'double');
        SHAPE.Zmin=fread(shpL,1,'double');
        SHAPE.Zmax=fread(shpL,1,'double');
        SHAPE.Mmin=fread(shpL,1,'double');
        SHAPE.Mmax=fread(shpL,1,'double');

    %% saving file pointers

    fseek(shxB,0,'eof');

    SHAPE.shpL = shpL;
    SHAPE.shxB = shxB;
    
    
    %% DBF data
    SHAPE.attributes = dbfread(FName,verbose)';
    
    % deblank fieldname (attribute name)
    for ia = 1:numel(SHAPE.attributes)
        SHAPE.attributes(ia).fieldname = deblank(SHAPE.attributes(ia).fieldname);
    end

    if verbose
        fprintf('\nAttributes of shape file: %s\n',FName);

        flds = fieldnames(SHAPE.attributes);
        flds = [flds(3:end); flds(2); flds(1)];
        fl = cellfun(@length,flds);
        fw = 4;
        s =  char(' '*ones(max(fl)+1,numel(flds)*fw));
        for ifld=1:numel(flds)
            s(max(fl)-fl(ifld)+1:max(fl),ifld*fw)=flds{ifld};
        end
        s(end,:)='=';     s = [s(end,:); s];
        for j=1:size(s,1)
            fprintf('%s\n',s(j,:));
        end

        fmtC = sprintf(' %%%ds',fw-1); % sprintf(' %%%s',fw);
        fmtN = sprintf(' %%%dd',fw-1); % sprintf(' %%%d',fw);

        NF = numel(flds);
        for ia = 1:numel(SHAPE.attributes)
            for ifld = 1:NF-3
                    fprintf(fmtN,SHAPE.attributes(ia).(flds{ifld}));
            end
            ifld = NF-2;
            fprintf(fmtN,numel(SHAPE.attributes(ia).(flds{ifld})));
            for ifld =NF-1:NF
                fprintf(fmtC,SHAPE.attributes(ia).(flds{ifld}));
            end
            fprintf('\n');
        end
    end

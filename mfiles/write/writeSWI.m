function writeSWI(basename,swi)
%WRITESWI writes input for SWI package (seawater intrusion package)
%
% Example:
%    writeSWI(basename,swi)
%
%    0.- Opens a new empty file to write all the data in
%    1.- Types in the new file the name and the date, does not work with this file
%    2.- Types in the command window the file which is being written
%    3.- Types the first line with the four requires values (NSRF ISTRAT ISWIZT
%        NPRN)
%    4.- Types the second line with the four required values (TOESLOPE TIPSLOPE
%        ZETAMIN DELZETA)
%    5.- Types the density values depending of the value given to ISTRAT
%    6.- The interfaces are written by number of plane and, afterwards, by layer
%    7.- Types the value of the effective porosity by layer
%    8.- Types the kind of source at each point by layer
%
%    The .swi file is ready to be used
%
% TO 090624

fid=fopen([basename,'.',swi.ext],'wt');                                     
%% 0

%fprintf(fid,'# MATLAB writeSWI %s\n',datestr(now));  % SWI doesn't allow comment lines
fprintf(    '# MATLAB writeSWI %s\n',datestr(now));                         

%% 1

fprintf(fid,'%10d%10d%10d%10d     NSRF ISTRAT ISWIZT NPRN\n',...            
        swi.NSRF,swi.ISTRAT,swi.ISWIZT,swi.NPRN);
    % number of inerfaces, swithc for stratified or interpolated and number
    % of steps between zeta recordings
    
%% 2 

fprintf(fid,'%10f%10f%10f%10f     TOESLOPE TIPSLOPE ZETAMIN DELZETA\n',...  
        swi.TOESLOPE,swi.TIPSLOPE,swi.ZETAMIN,swi.DELZETA);
        % max slop of toecell, of tipcell, minimum zeta before zeta is
        % removed from a cell, delzeta is elevation when it is moved into
        % an adjacent cell
%% 3
% Values for dimensionless density (NSRF+1 if ISTRAT=1, NSRF if ISTRT=0
warray(fid,swi.NU(:)',swi.unit,sprintf('(%dF15.5)',length(swi.NU)),'NU',true,swi.FREE);  

% the name of the ZETA field may be 'values' or 'term'
if isfield(swi.ZETA,'values')
    fieldNm='values';
else
    fieldNm='term';
end
    
if isstruct(swi.ZETA)
    % ZETA -- SWI wants for each surface and then for each layer the ZETA
    % this is compatible with the way data is stored in the budget file
    % in mflab. So if you have the struct with zeta planes in this way,
    % e.g obtained from reading a ZTA file:
    % ZETA=readbud([basename '.ZTA']);
    % then you could provide ZETA(it) as input to SWI in mfLab. it=end is default to
    % easily allow continuation from a previous run;
    
    swi.ZETA = swi.ZETA(1); % choose first element (always)

    % the field name of the arrays, one for each interface is 'values' or
    % 'term'. The labels are arbitrary. SWI uses ZETAPLANE1, ZETAPLANE2
    % etc. But these are not necessary for us as 'term' will be generated
    % by readBud when reading file [basename '.ZTA]
    
    % swi.NSRF must match the number of arrays in the struct
    if swi.NSRF ~= numel(swi.ZETA.(fieldNm))
        error('number of zeta planes =%d\ndoes not match the number of arrays in the struct swi.ZETA=%d',...
            swi.NSRF,numel(swi.ZETA.(fieldNm)));
    end
    
    if ~iscell(swi.ZETA.(fieldNm))
        error('swi.ZETA.%s must be a cell array, but is a %s',fieldNm,class(swi.ZETA.(fieldNm)));
    end
    
    % writing the zeta to the zetafile
    for ipln=1:swi.NSRF
        for iLay=1:swi.GRID.Nlay % for each layer
            Nlay = size(swi.ZETA.(fieldNm){ipln},3);
            z = swi.ZETA.(fieldNm){ipln}(:,:,min(iLay,Nlay));
            z = max(swi.GRID.ZBlay(:,:,iLay),min(swi.GRID.ZTlay(:,:,iLay),z));
            warray(fid,z,...
                swi.unit,'(10E15.6)',sprintf('Layer(%d), ZETA(%d)',iLay,ipln),true,swi.FREE);
        end
    end

elseif isnumeric(swi.ZETA)    
    % ZETA might be specified alternatively as a single array of size
    % NROW,NCOL,NSRF. Each layer in this array is the elevation of the ZETAPLANE
    % in 3D space.
    for ipln=1:swi.NSRF  % for each surface
        Nlay = size(swi.ZETA,3);
        for iLay=1:swi.GRID.Nlay % for each layer
            z = swi.ZETA(:,:,min(iLay,Nlay));
            z = max(swi.GRID.ZBlay(:,:,iLay),min(swi.GRID.ZTlay(:,:,iLay),z));
            warray(fid,z,...
                swi.unit,'(10E15.6)',sprintf('Layer(%d), ZETA(%d)',iLay,ipln),true,swi.FREE);
        end
    end
else
    error('mfLab:writeSWI:ZETA_wrong_structure',...
        ['ZETA should have same sturcture as if produced by readBud()\n',...
         'or be a Ny*Nx*NSRF array of doubles, each layer representing the\n',...
         'elevation of the ZETA plane in 3D space\n']);
end
%% 5
% Types the value of the effective porosity by layer
for iLay=1:swi.GRID.Nlay
    warray(fid,swi.SSZ(:,:,iLay)    ,swi.unit,'(12F12.3)',...
        sprintf('SSZ(%d)=effective porosity'    ,iLay),true,swi.FREE); 
end

%% 6
% Types the kind of source at each point by layer
for iLay=1:swi.GRID.Nlay
    warray(fid,swi.ISOURCE(:,:,iLay),swi.unit,'(25I4)',sprintf('ISOURCE(%d)',iLay),true,swi.FREE); 
end

%% The .swi file is ready to be used
fclose(fid);

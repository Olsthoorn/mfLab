function writeSWI2(basename,swi2)
%WRITESWI writes input for SWI package (seawater intrusion package)
%
% TO 090624

fid=fopen([basename,'.',swi2.ext],'wt');                                     
%% 0

%fprintf(fid,'# MATLAB writeSWI %s\n',datestr(now));  % SWI doesn't allow comment lines
fprintf(    '# MATLAB writeSWI2 %s\n',datestr(now));                         

%% 1

fprintf(fid,'# dataset 1\n');
fprintf(fid,'%10d%10d%10d%10d%10d%10d',...     % NSRF ISTRAT NOBS ISWIZT ISWIBD ISWIOBS',...            
        swi2.NSRF,  ...
        swi2.ISTRAT,...
        swi2.NOBS,  ...
        swi2.ISWIZT,...
        swi2.ISWIBD,...
        swi2.ISWIOBS);
    
if swi2.ADAPTIVE~=0
    fprintf(fid,'   %s','ADAPTIVE'); 
end

fprintf(fid,'\n');

%% 2a
fprintf(fid,'# dataset 2a\n');
fprintf(fid,'%10d%10d%10d     SOLVER, IPRSOL, MUTSOL\n',...  
        swi2.NSOLVER,  ...
        swi2.IPRSOL, ...
        swi2.MUTSOL);
        % max slop of toecell, of tipcell, minimum zeta before zeta is
        % removed from a cell, delzeta is elevation when it is moved into
        % an adjacent cell

if swi2.NSOLVER >1

    %% 2b
    fprintf(fid,'# dataset 2b\n');
    fprintf(fid,'%10d%10d%10d %12g %12g %12g %10d %12g',...
        swi2.MXITER, ...
        swi2.ITER1,  ...
        swi2.NPCOND, ...
        swi2.ZCLOSE, ...
        swi2.RCLOSE, ...
        swi2.RELAX,  ...
        swi2.NBPOL,  ...
        swi2.DAMP);

    if isfield(swi2,'DAMPT')
        fprintf(fid,'%12g',swi2.DAMPT);
    end
    fprintf(fid,'     MXITE ITER1 NPCOND ZCLOSE RCLOSE RELAC NBPOL DAMP [DAMPT]\n');
end

%% 3a 
fprintf(fid,'# dataset 3a\n');
fprintf(fid,'%12g %12g %12g %12g     TOESLOPE TIPSLOPE ALPHA BETA\n',...  
        swi2.TOESLOPE,swi2.TIPSLOPE,swi2.ALPHA,swi2.BETA);
        % max slop of toecell, of tipcell, minimum zeta before zeta is
        % removed from a cell, delzeta is elevation when it is moved into
        % an adjacent cell

%% 3b
if swi2.ADAPTIVE~=0
    fprintf(fid,'# dataset3b\n');
    fprintf(fid,'%10f%10f %12g\n',swi2.NADPTMX, swi2.NADPTMN, swi2.ADPTFCT);
end

%% 4 Values for dimensionless density (NSRF+1 if ISTRAT=1, NSRF if ISTRT=0

if swi2.ISTRAT<1, NNU = swi2.NSRF+2; else NNU = swi2.NSRF+1; end

fprintf(fid,'# dataset 4\n');
warray(fid,swi2.NU(1:NNU),swi2.unit,sprintf('(%dF15.5)',NNU),'',true,swi2.FREE);  

%% 5 ZETA

if isfield(swi2.ZETA,'values')
    fieldNm='values';
else
    fieldNm='term';
end

    
if isstruct(swi2.ZETA)
    % ZETA -- SWI wants for each surface(PLNE) and then for each layer the ZETA
    % this is compatible with the way data is stored in budget file
    % in mflab. So if you have the stuct with zeta planes in this way
    % ZETA=readbud([basename '.ZTA']);
    % ZETA=ZETA(i);
    % you could provide ZETA(i) as input to SWI in mfLab. i=end is default to
    % easily allow continuation from a previous run;

    swi2.ZETA = swi2.ZETA(1);
    
    if ~iscell(swi2.ZETA.(fieldNm))
        error('ZETA field %s must by cell array not %',fieldNm,class(swi2.ZETA.(fieldNm)));
    end
    
    for ipln=1:swi2.NSRF
        Nlay = size(swi2.ZETA.(fieldNm){ipln},3);
        for iLay=1:swi2.GRID.Nlay % for each layer
            z = swi2.ZETA.(fieldNm){ipln}(:,:,min(Nlay,iLay));
            z = max(swi2.GRID.ZBlay(:,:,iLay),min(swi2.GRID.ZTlay(:,:,iLay),z));
                %% 5
                warray(fid,z,...
                    swi2.unit,'(10E15.6)',sprintf('Layer(%d), ZETA(%d)',iLay,ipln),true,swi2.FREE);
        end
    end

elseif isnumeric(swi2.ZETA)    
    % ZETA might be specified alternatively as a single array of size
    % NROW,NCOL,NSRF. Each layer in this array is the elevation of the ZETAPLANE
    % in 3D space.
    for ipln=1:swi2.NSRF  % for each surface
        NLay = size(swi2.ZETA,3);
        for iLay=1:swi2.GRID.Nlay % for each layer
            z = swi2.ZETA(:,:,min(NLay,iLay));
            z = max(swi2.GRID.ZBlay(:,:,iLay),min(swi2.GRID.ZTlay(:,:,iLay),z));
            warray(fid,z,...
                swi2.unit,'(10E15.6)',sprintf('Layer(%d), ZETA(%d)',iLay,ipln),true,swi2.FREE);
        end
    end
else
    error('mfLab:writeSWI:ZETA_wrong_structure',...
        ['ZETA should have same sturcture as if produced by readBud()\n',...
         'or be a Ny*Nx*NSRF array of doubles, each layer representing the\n',...
         'elevation of the ZETA plane in 3D space\n']);
end

%% 6 SSZ (porosity)

for iLay=1:swi2.GRID.Nlay
    warray(fid,swi2.SSZ(:,:,iLay)    ,swi2.unit,'(12F12.3)',...
        sprintf('SSZ(%d)=effective porosity'    ,iLay),true,swi2.FREE); 
end

%% 7 ISOURCE

for iLay=1:swi2.GRID.Nlay
    warray(fid,swi2.ISOURCE(:,:,iLay),swi2.unit,'(25I4)',sprintf('ISOURCE(%d)',iLay),true,swi2.FREE); 
end

%% 8  OBSERVATIONS

% swi2.NOBS always zero in mfLab. It is a useless option.
if swi2.NOBS>0
    for iob=1:swi2.NOBS
        fprintf(fid,'%20s%10d%10d%10d\n',swi2.OBS(iob).name,siw2.OBS(iob).LRC);
    end
end
    

%% The .swi2 file is ready to be used
fclose(fid);

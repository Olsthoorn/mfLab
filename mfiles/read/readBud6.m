function B = readBud6(fname)
%READBUD6 reads compact budget file needed by MODPATH6
%
%  MODPATH release: Version 4.00 (V4, Release 1, 2-2000)
%    New routines to read MODFLOW-2000 budget files
%
% This fiele was translated from the original FORTRAN source code,
% with some simplifications and adpatations
%
% This function does not yet allow subselections of labels, stress periods,
% time steps etc.
%
% TO 130205
    
    fid = fopen(fname,'r');
    if fid<0
        error('%s: Can''t open budget file <<%s>>',mfilename,fname);
    end
    
    % find number of shorts to skip
    for nskip=0:4
        frewind(fid);
        fread(fid,nskip,'short');
        fread(fid,2,'int');
        lbl = fread(fid,16,'*char')';
        if all(ismember(lbl,[' ','A':'Z','0':'9','-+_'])),
            break;
        end
    end
    
    %nskip    = 2; % skip too shorts (record length indicator)
    
    [nRecOut,NMAXLBL] = readTheFile(fid,nskip); %#ok
    
    % fprintf('%s, nRecOut = %d, NMAXLBL = %d\n',mfilename,nRecOut,NMAXLBL);
    
    B = repmat(...
        struct(...
            'period' ,0,...
            'tstp'   ,0,...
            'label'  ,{[]},...
            'NLAY'   ,0,'NROW',0,'NCOL',0,...
            'term'   ,{[]},...
            'lays'   ,0,...
            'rows'   ,0,...
            'cols'   ,0,...
            'pertim' ,0,...
            'totim'  ,0,...
            'delt'   ,0,...
            'NIFACE' ,0,...
            'CTMP',  {[]}),...
        nRecOut, 1);

    % preallocate
    B = readTheFile(fid,nskip,B);
    
    fclose(fid);
end

function varargout = readTheFile(fid,nskip,B)
    % DIMENSION VAL(20)
    % DIMENSION B(nRecOut).term{iLbl}(NCOL,NROW,NLAY),IBUFF(NCOL,NROW,NLAY)
    % CHARACTER*16 LABEL

    %....THIS ROUTINE READS THROUGH BUDGET FILE AND RETURNS THE NEXT SET OF
    %....DATA THAT MATCHES THE STRESS PERIOD AND TIME STEP. THIS ROUTINE READS
    %....DATA INTO B(nRecOut).term{iLbl} EVEN if THE DATA ARE POINT DATA.
    %....Budget file types as indicated in NBTYPE
    %      0 for 3-D array in original format
    %      1 for 3-D array with extra header
    %      2 for list without auxiliary variables
    %      3 for 1-layer data array and 1-layer indicator array
    %      4 for layer-1 array
    %      5 for list with auxiliary data
    
    verbose   = false;
    SCANNING = nargin<3;
    
    if SCANNING, mode='Scanning'; else mode='Reading'; end %#ok
    
    nRecIn     = 0;
    nRecOut    = 0;
    KPold      = 0; % period of previous record
    KSold      = 0; % kstep  of previous record
    NMAXLBL        = 0; % highest label number on any record

    [precision,pbytes,NLAY,NROW,NCOL] = budgetPrecision(fid,nskip);
    
    frewind(fid);

    NCELLS = NROW*NCOL*NLAY;
    NRC    = NROW*NCOL;

    while ~feof(fid)
        nRecIn  = nRecIn +1;

        [KS,KP,LABEL,~,~,NLCODE] = firstHdr(fid,nskip);
        if feof(fid)
            break;
        end
  
        COMPACT = NLCODE<0;

        if KP>KPold || (KP==KPold && KS>KSold)
            nRecOut = nRecOut+1;
            iLbl    = 1;
            if ~SCANNING
                B(nRecOut).period = KP;
                B(nRecOut).tstp   = KS;
                B(nRecOut).NCOL   = NCOL;
                B(nRecOut).NROW   = NROW;
                B(nRecOut).NLAY   = NLAY;
                B(nRecOut).label{iLbl} = LABEL;
                B(nRecOut).term{iLbl}  = zeros(NCOL,NROW,NLAY);
            end
        else
            iLbl    = iLbl + 1;
            if SCANNING
                NMAXLBL = max(NMAXLBL,iLbl);
            else
                B(nRecOut).label{iLbl} = LABEL;
            end
        end

        % Specific for compact budget file
        if COMPACT
            if SCANNING
                NBTYPE = secondHdr(fid,nskip,precision);
            else
            [NBTYPE,...
             B(nRecOut).delt,...
             B(nRecOut).pertim,...
             B(nRecOut).totim] = secondHdr(fid,nskip,precision);
            end

            if NBTYPE==5 
                fread(fid,nskip,'short');
                NVAL = fread(fid,1,'int');
                fread(fid,nskip,'short');
                %IBYTES=IBYTES + 4
                if NVAL>1
                    B(nRecOut).CTMP{NVAL-1} = 'dummy';
                    fread(fid,nskip,'short');
                    for iv = 1:NVAL-1
                        B(nRecOut).CTMP{iv} = fread(fid,16,'*char');
                        % IBYTES=IBYTES + 16*(NVAL-1)
                    end
                    fread(fid,nskip,'short');
                    for iv=1:NVAL-1
                        if strcmp(B(nRecOut).CTMP{iv},'IFACE'),
                            B(nRecOut).NIFACE=iv+1;
                        end
                    end
                end
            end
            
            switch NBTYPE
                case {0,1}
                    fread(fid,nskip,'short');
                    if SCANNING
                       fseek(fid,pbytes*NCELLS,0);
                    else
                        B(nRecOut).term{iLbl}(:,:,:) = ...
                            reshape(fread(fid,NCELLS,precision),[NCOL,NROW,NLAY]);
                    end
                    fread(fid,nskip,'short');
                case {2,5}
                    fread(fid,nskip,'short');
                    NLST = fread(fid,1,'int');
                    fread(fid,nskip,'short');
                    B(nRecOut).term{iLbl} = zeros(NCOL,NROW,NLAY);
                    if NLST>0
                        for N=1:NLST
                            fread(fid,nskip,'short');
                            if SCANNING
                                fseek(fid, 4    ,0);
                                fseek(fid,pbytes,0);
                            else
                                ICELL   = fread(fid,1,'int');
                                VAL     = fread(fid,1,precision);
                                k       = floor( (ICELL-1)/NRC) + 1;
                                i       = floor(((ICELL - (k-1)*NRC)-1 )/NCOL) + 1;
                                j       = ICELL - (k-1)*NRC - (i-1)*NCOL;
                                B(nRecOut).term{iLbl}(j,i,k)=B(nRecOut).term{iLbl}(j,i,k)+VAL(1);
                            end
                            fread(fid,nskip,'short');
                        end
                    end
                case 4
                    fread(fid,nskip,'short');
                    if SCANNING
                        fseek(fid,pbytes*NRC,0);
                    else
                        B(nRecOut).term{iLbl}        = zeros(NCOL,NROW,NLAY);
                        B(nRecOut).term{iLbl}(:,:,1) = fread(fid,[NCOL,NROW],precision);
                    end
                    fread(fid,nskip,'short');
                case 3
                    fread(fid,nskip,'short');
                    if SCANNING
                        fseek(fid,     4*NRC,0);                        
                        fread(fid,nskip,'short');
                        fread(fid,nskip,'short');
                        fseek(fid,pbytes*NRC,0);
                    else
                        [Ix,Iy] = meshgrid(1:NCOL,(1:NROW)');
                        Iz      = fread(fid,[NCOL,NROW],'int');
                        Idx = cellIndex(Ix(:),Iy(:),Iz(:),[NROW,NCOL,NLAY]);
                        fread(fid,nskip,'short');
                        fread(fid,nskip,'short');
                        V      = fread(fid,[NCOL,NROW],precision);
                        B(nRecOut).term{iLbl}=zeros(NCOL,NROW,NLAY);
                        B(nRecOut).term{iLbl}(Idx) = V(:);
                    end
                    fread(fid,nskip,'short');
                otherwise
                    error('%s: Unknown value NBTYPE = %d, must be between 0 and 5',mfilename,NBTYPE);
            end

        else
            fread(fid,nskip,'short');
            if SCANNING
                fseek(fid,pbytes*NCELLS,0);
            else
                B(nRecOut).term{iLbl} = ...
                    reshape(fread(fid,NCELLS,precision),[NCOL,NROW,NLAY]);
            end
            fread(fid,nskip,'short');
        end

        if ~SCANNING % convert to MATLAB native ordering
            B(nRecOut).term{iLbl} = permute(B(nRecOut).term{iLbl},[2,1,3]);
        end
        KPold = KP;
        KSold = KS;
    
        % print info for the user who's waiting
        if verbose
            fprintf(2,'%s %s for Period %d, Time Step %d\n',mode,LABEL,KP,KS); %#ok
        else
            fprintf('.');
            if rem(nRecIn,50)==0
                fprintf('%d\n',nRecOut)
            end
        end
        
    end

    
    fprintf('%d\n',nRecIn);

    if SCANNING
        fprintf('%s: %d input records scanned,  %d output records in budget file\n',...
            mfilename,nRecIn,nRecOut);
        varargout{1} = nRecOut;
        varargout{2} = NMAXLBL;
    else
        fprintf('%s: %d input records read yielding  %d elements in output struct\n',...
        mfilename,nRecIn,nRecOut);
        varargout{1} = B;
    end
end
     
function [KSTP,KPER,LABEL,NC,NR,NL] = firstHdr(fid,nskip)    
    fread(fid,nskip,'short');
    KSTP  = fread(fid,1,'int');
    KPER  = fread(fid,1,'int');
    LABEL = fread(fid,16,'*char')';
    NC    = fread(fid,1,'int');
    NR    = fread(fid,1,'int');
    NL    = fread(fid,1,'int');
    fread(fid,nskip,'short');
end

function [ICODE,DELT,PERTIM,TOTIM] = secondHdr(fid,nskip,precision)
    fread(fid,nskip,'short');
    ICODE   = fread(fid,1,'int');
    DELT    = fread(fid,1,precision);
    PERTIM  = fread(fid,1,precision);
    TOTIM   = fread(fid,1,precision);
    fread(fid,nskip,'short');
end

function [precision,pbytes,NLAY,NROW,NCOL] = budgetPrecision(fid,nskip)
    %  Determine single or double precision file type for a MODFLOW
    %  budget file:  0=unrecognized, 1=single, 2=double.
    % USE DOUBLEBUDGET, ONLY:   IPRBUD,DBLBUFF,SNGBUFF
    % DIMENSION B(nRecOut).term{iLbl}(NCOL,NROW,NLAY)
    % DOUBLE PRECISION DELTD,PERTIMD,TOTIMD,VALD
    % REAL(KIND=4) DELTS,PERTIMS,TOTIMS,VALS
    % CHARACTER*16 TEXT1,TEXT2
    % INTEGER IOUT
    
    frewind(fid);
    
    precisions = {'single','double','undefined'};

    [~,~,LABEL1,NCOL,NROW,NLCODE] = firstHdr(fid,nskip);

    NLAY = abs(NLCODE);
    
    COMPACT = NLCODE<0;
    
    NCELLS = NCOL*NROW*NLAY;
    
    fp = ftell(fid); % remember filepos
    
    for ip = 1:3
        precision = precisions{ip};
        %  Read data depending on ICODE.  ICODE 0,1, or 2 are the only allowed
        %  values because the first budget terms must be from the internal
        %  flow package (BCF,LPF, or HUF).

        fseek(fid,fp,-1); % reset file pointer
        
        if COMPACT
            ftype ='Compact';
            ICODE = secondHdr(fid,nskip,precision);
            switch ICODE
                case {0,1}
                    fread(fid,nskip,'short');
                    fread(fid,NCELLS,precision);
                    fread(fid,nskip,'short');
                case 2
                    fread(fid,nskip,'short');
                    NLST = fread(fid,1,'int');
                    fread(fid,nskip,'short');

                    if NLST>0
                        for N=1:NLST
                            fread(fid,nskip,'short');
                            fread(fid,1,'int'); % ICELL
                            fread(fid,1,precision); % VAL
                            fread(fid,nskip,'short');
                        end
                    end
                otherwise
                    error('%s: Illegal ICODE = %d, must be 0,1 or 2',mfilename,ICODE);
            end
        else
            ftype = 'Full';
            % ICODE == 0 (old full budget file)
            fread(fid,nskip,'short');
            fread(fid,NCELLS,precision);
            fread(fid,nskip,'short');
        end

        %  Read 2nd header and check for valid type.
        [~,~,LABEL2] = firstHdr(fid,nskip);

        legalChars = [' ','A':'Z','0':'9','-+_()'];
        
        if all(ismember(LABEL1,legalChars)) && all(ismember(LABEL2,legalChars)) 
%         if (strcmpi(LABEL1,'         STORAGE') && strcmpi(LABEL2,'   CONSTANT HEAD')) || ...
%            (strcmpi(LABEL1,'   CONSTANT HEAD') && strcmpi(LABEL2,'FLOW RIGHT FACE '))
            switch precision
                case 'single', break;
                case 'double', break;
            end
        end
    end

    fprintf('\n');

    switch precision
        case 'single'
            fprintf(1,'This is a Single Precision Binary %s Budget File\n',ftype);
            pbytes = 4;
        case 'double'
            fprintf(1,'This is a Double Precision Binary  %s Budget file\n',ftype);
            pbytes = 8;
        otherwise
            error('%s: Unable to determine the precision of the budget file',mfilename);
    end
end

function nam = writeNAM(basename,nam)
%WRITENAM writes name files for MODFLOW, MT3DMS and SEAWAT
%
% Example:
%    writeNAM(basename,namfiles)
%
% TO 070630 100120 130203

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

    fmt='%-15s%5d %s.%s\n';         % for nam file
    if ispc
        fms='"%s\\%s%s" %s.nam\n';  % for batfile
    else
        fms='"%s/%s%s" %s.nam\n';   % for batfile
    end

    for iModel=1:length(nam.MODEL)                
        
        % Put LST in front
        i = strmatchi('LIST',nam.PCKG,'exact');
        nam.EXT{i} = [nam.MODEL{iModel} '.LST'];
        if any(i)
            i = i(1);
            I = 1:numel(nam.PCKG);
            I = I([i,1:i-1,i+1:end]);
            nam.PCKG = nam.PCKG(I);
            nam.UNIT = nam.UNIT(I);
            nam.EXT  = nam.EXT(I);
            nam.SCEN = nam.SCEN(I);
        end

        hdr=sprintf('# MATLAB       %s writeNAM --> nam file for model %s\n', datestr(now),nam.MODEL{iModel});
            fprintf('# MODFLOW etc. %s writeNAM --> nam file for model %s\n', datestr(now),nam.MODEL{iModel});

        % get file base name, extension and path of the models for which to generate input 

        [MP, MF , ME]=fileparts(nam.mdlpath{iModel});
        fid =fopen([MF, '.nam'],'wt');
        fprintf(fid ,'%s',hdr);

        switch nam.MODEL{iModel}
            % 1:  MF bat and nam file
            case 'MODPATH'  % version 6 only
                
                % add MPBAS to the PCKG and the nam struct to return it to mf_setup
                % for later use with writeMPBAS
                nam.UNIT = [nam.UNIT; nextFreeUnits(nam.UNIT,1)];
                nam.EXT  = [nam.EXT ; 'MPBAS'];
                nam.PCKG = [nam.PCKG; 'MPBAS'];

                fprintf(fid,fmt,'MPBAS',nam.UNIT(end),basename,nam.EXT{end});
                
                I = strmatchi('DIS',nam.PCKG);
                if ~I(1)
                    error('%s: no DIS file, switch on ''DIS'' in the nam file',mfilename);
                else
                    fprintf(fid,fmt,nam.PCKG{I(1)},nam.UNIT(I(1)),basename,nam.EXT{I(1)});
                end
                
                % head and budget file
                I = strmatchi('DATA(BINARY)',nam.PCKG);
                for i=I(:)'
                    switch nam.EXT{i}
                        case 'HDS'
                            fprintf(fid,fmt,'HEAD',nam.UNIT(i),basename,nam.EXT{i});
                        case {'BGT','BUD'}
                            fprintf(fid,fmt,'BUDGET',nam.UNIT(i),basename,nam.EXT{i});
                        otherwise
                    end
                end
                
            otherwise % i.e. MODFLOW, MT3DMS, SEAWAT etc but not MODPATH
                % Use PCKG in nam.PCKG verbatim
                for i=1:length(nam.PCKG)
                    if strmatchi(nam.PCKG{i},nam.LegalPCKG{iModel})
                        fprintf(fid,fmt,nam.PCKG{i},nam.UNIT(i),basename,nam.EXT{i});
                        fprintf(    fmt,nam.PCKG{i},nam.UNIT(i),basename,nam.EXT{i});
                    end
                end
        end
        fclose(fid);

        % mf batfile
        fid=fopen([MF,'.bat'],'wt');
        fprintf(fid,fms,MP,MF,ME,MF);
        fclose(fid);
    end
end

function units = nextFreeUnits(UNIT,n)
    %% unit= nextFreeUnits(UNIT,n) -- find next free unit numbers (Below 89)

    UNIT = sort(UNIT);          % occupied unit numbers
    
    units = (10:100)';          % start with unit 10 to prevent STDIN (5) and STDOUT (6)
    
    free = ~ismember(units,UNIT);
    
    units = units(free);
    units = units(1:n);            
end
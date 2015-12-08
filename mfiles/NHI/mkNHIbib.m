function mkNHIbib(NHIzip,NHIbib)
%MKNHIBIB starts with directory with NHIzipfiles and generates a directory with unzipped files
%
% Example:
%    mkNHIbib(NHIzip,NHIbib)
%
% Starts with directory with NHIzipfiles and generate a directory
% with the unzipped files each set in its own subdir, skipping
% if subdir already exists.
%
% to renew remove directory with specific fileset
%
% TO 090727


%NHIbib='..\Bibliotheek_NHI\';
%AGVbib='..\Bibliotheek_AGV\';
%NHIzip='Z:\tolsthoorn On My Mac\GRWMODELS\NHI_and_AGV\ZIP090727\';

here=pwd;

cd(   NHIzip);
d=dir(NHIzip);

fprintf('\n\n\nzipName={...\n');
for i=1:length(d)
    if length(d(i).name)>4 && strcmp(d(i).name(end-3:end),'.zip')
        Newdir=[NHIbib d(i).name(1:end-4)];
        if ~exist(Newdir,'dir')
            fprintf('mkdir ''%s''\n',Newdir);
            fprintf('unzip ''%s'' in ''%s''\n',d(i).name,Newdir);
            fnames=unzip(d(i).name,Newdir);
        else
            fprintf('skpping aready existing dir ''%s''\n',Newdir);
        end
    end
end

cd(here)

%% 
     

%
% 
% bibName={...
% 'anisotropie_factor';...
% 'anisotropie_hoek';...
% 'beregeningslokaties';...
% 'bodemfysische_eenheid';...
% 'bodemhoogte_hoofd_';...
% 'bodemhoogte_prim_';...
% 'bodemhoogte_sec_';...
% 'bodemhoogte_tert_';...
% 'bot_laag';...
% 'buisdrainage_bodh';...
% 'buisdrainage_c';...
% 'c_laag';...
% 	};



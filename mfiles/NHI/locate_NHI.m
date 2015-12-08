%LOCATE_NHI starts with directory with NHIzipfiles and generate a directory and ...
% with the unzipped files each set in its own subdir
% skips if subdir already exists.
% to renew remove directory with specific fileset
%
% TO 090727

AGVbib='..\Bibliotheek_AGV\';
NHIzip='Z:\tolsthoorn On My Mac\GRWMODELS\NHI_and_AGV\ZIP090727\';
NHIbib='Z:\tolsthoorn On My Mac\GRWMODELS\NHI_and_AGV\Bibliotheek_NHI\';

%% Unzip the NHI model into NHIbib
mkNHIbib(NHIzip,NHIbib)

%% 
% NHI2AGV
% extracts AGV model from NHI datafiles
% JB 090722

%base='D:\data\NHI\model\input\mf\v067\';
base='Z:\tolsthoorn On My Mac\GRWMODELS\MYWORK\AGVmodel\';
home=[base 'AGV1'];

path(path,home);

cd(home);
fprintf('Running NHI2AGV from dir ''%s''\n',home);
fprintf('%s\n',datestr(now));

%% Initialize size of AGV model
xLimAGV=[100000 155000];
yLimAGV=[450000 500000];

%xLimAGV=[110000 145000];
%yLimAGV=[460000 490000];

%xLimAGV=[120000 135000];
%yLimAGV=[460000 480000];

%xLimAGV=[132500 135000];
%yLimAGV=[465000 467500];

% =================== READ ASC FILES ===================== 
fprintf('Determining size of NHI model from kDlaag director in NHIbib\n');
fprintf('''%s''\n',NHIbib);

%% Determine size of AGV model

d=dir(NHIbib);  % ascii NHI data set files
for i=1:length(d)
    if  strcmp(d(i).name,'dikte_laag');  % look for kd_laag directory in bib
        dsub=dir([NHIbib d(i).name '\dikte_laag*.asc']);
        NLAY=length(dsub);
        FName=[NHIbib d(i).name '\' dsub(end).name]; % this is dikte_laag(end).asc
        [dummy,xGr,yGr,xm,ym,IxLim,IyLim]=readASC(FName,'',xLimAGV,yLimAGV);
        NCOL=length(xm);
        NROW=length(ym);
        fprintf('Number of layers according to last "dikte_laag" = %c\n',...
            dsub(end).name(end-4));
        break; % ready
    end
end  

fprintf('AGV model: NCOL,NROW,NLAY= %d %d %d\n',NCOL,NROW,NLAY);


%% Allocate space for matrices



%% Read the file, deduce filetype from meta data
d=dir(NHIbib);

tic

for i=1:length(d)
    if d(i).name(1)~='.' && d(i).isdir
        fprintf('%s\n',d(i).name);
        dsub=dir([NHIbib d(i).name '\*.met']);
        FNmeta=[NHIbib d(i).name '\' dsub(1).name];
        meta=readNHImeta(FNmeta);

        filetype={'ArcInfo_shapefile','ASCII_MODFLOW','ArcInfo_ASCII_raster'};

        switch meta.Type
            case filetype(1)
%                 pth=[NHIbib d(i).name];
%                 setshapedir([pth,'\'],d(i).name)
%                 dsub=dir([NHIbib d(i).name '\*.SHP']);
%                 switch meta.Zip_file(1:end-4)
%                     case {'districten'}  % skip or draw
%                         districten=readshp([pth '\' dsub(1).name]);
%                         figure;
%                         plotshp(districten);
%                     case {'kenmerken_lsws'} % skip or draw
%                         kenmerken_lsws=readshp([pth '\' dsub(1).name]);
%                         figure;
%                         plotshp(kenmerken_lsws);
%                     case {'ligging_dm'} % skip or draw
%                        ligging_dm=readshp([pth '\' dsub(1).name]);
%                         figure;
%                         plotshp(ligging_dm);
%                  end
            case filetype(2)
                if ~exist('WEL','var')
                    pth=[NHIbib d(i).name];               
                    dsub=dir([NHIbib d(i).name '\*.SCD']);
                    switch meta.Zip_file(1:end-4)
                        case {'onttrekkingen'}
                            WEL=readSCD(dsub.name,pth,'w',IxLim,IyLim);
                    end
                end
            case filetype(3)
                pth=[NHIbib d(i).name];
                dsub=dir([NHIbib d(i).name '\*.ASC']);
                basename=meta.File_prefix;
                switch basename
                    case {'anisotropie_factor'}  % layer 1-3
                        if ~exist(basename,'var')
                            anisotropie_factor=readASC(dsub.name,pth,xLimAGV,yLimAGV);
                        end
                    case {'anisotropie_hoek'}    % layer 1-3
                        if ~exist(basename,'var')
                            anisotropie_hoek  =readASC(dsub.name,pth,xLimAGV,yLimAGV);
                        end
%                     case {'beregeningslokaties'}
%                         beregeningslokaties=readASC(dsub.name,pth,xLimAGV,yLimAGV);
%                     case {'bodemfysische_eenheid'}
%                         bodemfysische_eenheid=readASC(dsub.name,pth,xLimAGV,yLimAGV);
                    case {'bodemhoogte_hoofd_'}
                        if ~exist([basename 'j'],'var')
                            for j=1:length(dsub)
                                if strcmp(dsub(j).name,'_j')
                                    bodemhoogte_hoofd_j=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
                                end
    %                             if strcmp(dsub(j).name,'_w')
    %                                 bodemhoogte_hoofd_w=readASC(d(j).name,pth,xLimAGV,yLimAGV);
    %                             end
    %                             if strcmp(dsub(j).name,'z');
    %                                 bodemhoogte_hoofd_z=readASC(d(j).name,pth,xLimAGV,yLimAGV);
    %                             end
                            end
                        end
                    case {'bodemhoogte_prim_'}
                        if ~exist([basename 'j'],'var')
                            for j=1:length(dsub)
                                if ~isempty(findstr('_j',dsub(j).name))
                                    bodemhoogte_prim_j=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
                                end
    %                             if ~isempty(findstr('_w',dsub(j).name))
    %                                 bodemhoogte_prim_w=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
    %                             end
    %                             if ~isempty(findstr('_z',dsub(j).name))
    %                                 bodemhoogte_prim_z=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
    %                             end
                            end
                        end
%                     case {'bodemhoogte_sec_'}
%                         for j=1:length(dsub)
%                             if ~isempty(findstr('_j',dsub(j).name))
%                                 bodemhoogte_sec_j=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
%                             end
%                             if ~isempty(findstr('_w',dsub(j).name))
%                                 bodemhoogte_sec_w=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
%                             end
%                             if ~isempty(findstr('_z',dsub(j).name))
%                                 bodemhoogte_sec_z=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
%                             end
%                         end
%                     case {'bodemhoogte_tert_'}
%                         for j=1:length(dsub)
%                             if ~isempty(findstr('_j',dsub(j).name))
%                                 bodemhoogte_tert_j=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
%                             end
%                             if ~isempty(findstr('_w',dsub(j).name))
%                                 bodemhoogte_tert_w=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
%                             end
%                             if ~isempty(findstr('_z',dsub(j).name))
%                                 bodemhoogte_tert_z=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
%                             end
%                         end
                    case {'bot_laag'}  % bottom of layers in cm (is converted to m)
                        if ~exist('BOT','var')
                            BOT=NaN(NROW,NCOL,NLAY);
                            for iFile=1:length(dsub)
                                for iLay=1:NLAY
                                    Laynum=sprintf('%d',iLay);
                                    if ~isempty(findstr(Laynum,dsub(iFile).name))
                                        BOT(:,:,iLay)=readASC(dsub(iLay).name,pth,xLimAGV,yLimAGV)/100;
                                    end
                                end
                            end
                        end
                    case {'buisdrainage_bodh'}
                        if ~exist(basename,'var')
                            buisdrainage_bodh=readASC(dsub.name,pth,xLimAGV,yLimAGV);
                        end
                    case {'buisdrainage_c'}
                        if ~exist(basename,'var')
                            buisdrainage_c   =readASC(dsub.name,pth,xLimAGV,yLimAGV);
                        end
                    case {'c_laag'}
                        if ~exist('VCONT','var');
                            MINC=1; %minimum resistance of layers
                            VCONT =NaN(NROW,NCOL,NLAY-1);
                            for iFile=1:length(dsub)
                                for iLay=1:NLAY-1
                                    Laynum=sprintf('%d',iLay);
                                    if ~isempty(findstr(Laynum,dsub(iFile).name))
                                        VCONT(:,:,iLay)=1./max(MINC,readASC(dsub(iFile).name,pth,xLimAGV,yLimAGV));
                                    end
                                end
                            end
                        end
                    case {'cdr_hoofd'}
                        if ~exist(basename,'var')
                            cdr_hoofd=readASC(dsub.name,pth,xLimAGV,yLimAGV);
                        end
%                     case {'cdr_prim'}
%                         if isempty(dsub)
%                            fprintf('Skipping ''%s'', no data files found\n',d(i).name)
%                         else
%                            cdr_prim=readASC(dsub.name,pth,xLimAGV,yLimAGV);
%                         end
%                     case {'cdr_tert'}
%                         if isempty(dsub)
%                            fprintf('Skipping ''%s'', no data files found\n',d(i).name)
%                         else
%                             cdr_tert=readASC(dsub.name,pth,xLimAGV,yLimAGV);
%                         end
%                     case {'cinf_hoofd'}
%                         if isempty(dsub)
%                            fprintf('Skipping ''%s'', no data files found\n',d(i).name)
%                         else
%                             cinf_hoofd=readASC(dsub.name,pth,xLimAGV,yLimAGV);
%                         end
%                     case {'cinf_prim'}
%                         if isempty(dsub)
%                            fprintf('Skipping ''%s'', no data files found\n',d(i).name)
%                         else
%                             cinf_prim=readASC(dsub.name,pth,xLimAGV,yLimAGV);
%                         end
%                     case {'cinf_sec'}
%                         if isempty(dsub)
%                            fprintf('Skipping ''%s'', no data files found\n',d(i).name)
%                         else
%                             cinf_sec=readASC(dsub.name,pth,xLimAGV,yLimAGV);
%                         end
%                     case {'clconc_laag'}
%                         clconc=NaN(NROW,NCOL,NLAY);
%                         for iFile=1:length(dsub)
%                             for iLay=1:NLAY
%                                 Laynum=sprintf('%d',iLay);
%                                 if ~isempty(findstr(Laynum,dsub(iFile).name))
%                                     clconc(:,:,iLay)=readASC(dsub(iFile).name,pth,xLimAGV,yLimAGV);
%                                 end
%                             end
%                         end
                    case {'dikte_laag'}
                        if ~exist('THK','var')
                            THK=NaN(NROW,NCOL,NLAY);
                            for iFile=1:length(dsub)
                                for iLay=1:NLAY
                                    Laynum=sprintf('%d',iLay);
                                    if ~isempty(findstr(Laynum,dsub(iFile).name))
                                        THK(:,:,iLay)=readASC(dsub(iFile).name,pth,xLimAGV,yLimAGV)/100;
                                        % Divide by 100 as THK is given in cm
                                    end
                                end
                            end
                        end
                    case {'kd_laag'}
                        if ~exist('TRAN','var')
                            TRAN=NaN(NROW,NCOL,NLAY);
                            for iFile=1:length(dsub)
                                for iLay=1:NLAY
                                    Laynum=sprintf('%d',iLay);
                                    if ~isempty(findstr(Laynum,dsub(iFile).name))
                                        TRAN(:,:,iLay)=readASC(dsub(iFile).name,pth,xLimAGV,yLimAGV);
                                    end
                                end
                            end
                        end
                    case {'nl_001_'}  % Startheads
                        if ~exist('STRTHD','var')
                            STRTHD=NaN(NROW,NCOL,NLAY);
                            for iFile=1:length(dsub)
                                for iLay=1:NLAY
                                    Laynum=sprintf('%d',iLay);
                                    if ~isempty(findstr(['_0' Laynum],dsub(iFile).name))
                                        STRTHD(:,:,iLay)=readASC(dsub(iFile).name,pth,xLimAGV,yLimAGV);
                                    end
                                end
                            end
                        end
                    case {'oppervlakte_ow'}
                        if ~exist([basename 'h'],'var')
                            for j=1:length(dsub)
                                Var=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
                                if ~isempty(findstr('owh',dsub(j).name))
                                    oppervlakte_owh=Var;
                                end
                                if ~isempty(findstr('owp',dsub(j).name))
                                    oppervlakte_owp=Var;
                                end
    %                             if ~isempty(findstr('ows',dsub(j).name))
    %                                 oppervlakte_ows=Var;
    %                             end
    %                             if ~isempty(findstr('owt',dsub(j).name))
    %                                 oppervlakte_owt=Var;
    %                             end
                            end
                        end
                        case {'peil_hoofd_'}
                            if ~exist([basename 'j'],'var')
                                peil_hoofd=NaN(NROW,NCOL);
                                for j=1:length(dsub)
                                    if ~isempty(findstr('_j',dsub(j).name))
                                        peil_hoofd_j=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
                                    end
    %                                 if ~isempty(findstr('_w',dsub(j).name))
    %                                     peil_hoofd_w=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
    %                                 end
    %                                 if ~isempty(findstr('_z',dsub(j).name))
    %                                     peil_hoofd_z=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
    %                                 end
                               end
                            end
                         case {'peil_prim_'}
                             if ~exist([basename 'j'],'var')
                                for j=1:length(dsub)
                                    if ~isempty(findstr('_j',dsub(j).name))
                                        peil_prim_j=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
                                    end
    %                                 if ~isempty(findstr('_w',dsub(j).name))
    %                                     peil_prim_w=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
    %                                 end
    %                                 if ~isempty(findstr('_z',dsub(j).name))
    %                                     peil_prim_z=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
    %                                 end
                                end
                             end
%                         case {'peil_sec_'}
%                             for j=1:length(dsub)
%                                 if ~isempty(findstr('_j',dsub(j).name))
%                                     peil_sec_j=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
%                                 end
%                                 if ~isempty(findstr('_w',dsub(j).name))
%                                     peil_sec_w=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
%                                 end
%                                 if ~isempty(findstr('_z',dsub(j).name))
%                                     peil_sec_z=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
%                                 end
%                             end
%                         case {'peil_tert_'}
%                             for j=1:length(dsub)
%                                 if ~isempty(findstr('_j',dsub(j).name))
%                                     peil_tert_j=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
%                                 end
%                                 if ~isempty(findstr('_w',dsub(j).name))
%                                     peil_tert_w=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
%                                 end
%                                 if ~isempty(findstr('_z',dsub(j).name))
%                                     peil_tert_z=readASC(dsub(j).name,pth,xLimAGV,yLimAGV);
%                                 end
%                             end
                        case {'rch'}  % layer 1-3
                            if ~exist('RCH','var')
                                RCH(1).values=readASC(dsub.name,pth,xLimAGV,yLimAGV);
                            end
                        case {'top_laag'} % in cm is converted to m
                            if ~exist('TOP','var')
                                TOP=NaN(NROW,NCOL,NLAY);
                                for iFile=1:length(dsub)                              
                                    for iLay=1:NLAY
                                        Laynum=sprintf('%d',iLay);
                                        if ~isempty(findstr(Laynum,dsub(iFile).name))
                                            TOP(:,:,iLay)=readASC(dsub(iLay).name,pth,xLimAGV,yLimAGV)/100;
                                        end
                                    end
                                end
                            end
%                         case {'wortelzonedikte'}
%                             wortelzonedikte=readASC(dsub.name,pth,xLimAGV,yLimAGV);
                end
        end
    end
end

     
%%  Issue spent time
     toc
%% get types of files etc from NHVbib Practica tool
% filetype={'ASCII_MODFLOW','ArcInfo_ASCII_raster','ArcInfo_shapefile'};
% ift=2;
% for i=1:length(d)
%     if d(i).name(1)~='.' && d(i).isdir
%         dsub=dir([NHIbib d(i).name '\*.met']);
%         FNmeta=[NHIbib d(i).name '\' dsub(1).name];
%         meta=readNHImeta(FNmeta);
%         if strcmp(meta.Type,filetype{ift})
%             fprintf('case {''%s''}\n',meta.Zip_file(1:end-4));
%         end
%     end
% end 

%% Genenral head boundaries
% creating them from the ESRI asci rasters that have been read in

[DX,DY]=meshgrid(abs(diff(xGr)),abs(diff(yGr)));


%% Drain  -- using buisdrainage_bodh, buisdrainage_c
IDRN=find(~isnan(buisdrainage_bodh) & ~isnan(buisdrainage_c));
DRNindices=cellIndices(IDRN,[NROW,NCOL,NLAY]);
%plot(xm(DRNindices(:,2)),ym(DRNindices(:,1)),'.')
DRN=[DRNindices(:,[3,1,2]),buisdrainage_bodh(IDRN),DX(IDRN).*DY(IDRN)./buisdrainage_c(IDRN)];


%% GHB -- using peil_hoofd_j cdr_hoofd and oppervlakte_owh
% (bodemhoogte_hoofd was not among the data files in NHI website)

%peil_hoofd_j  % jaargemiddeld peil
%peil_hoofd_w  % winter peil 
%peil_hoofd_z  % zomer peil

IGHB=find(~isnan(cdr_hoofd) & ~isnan(peil_hoofd_j));
GHBindices=cellIndices(IGHB,[NROW,NCOL,NLAY]);
%plot(xm(GHBindices(:,2)),ym(GHBindices(:,1)),'.');
GHB=[GHBindices(:,[3,1,2]),peil_hoofd_j(IGHB),oppervlakte_owh(IGHB)./cdr_hoofd(IGHB)];

%% RIV using bodemhoogte_prim_j, peil_prim_j, oppervlakte_owp.
% cdr_prim was not in the data, we use a default of 100 days

cdr_prim_default=100; %days

%bodemhoogte_prim_j % jaar peil
%bodemhoogte_prim_w % winter peil
%bodemhoogte_prim_z % zomer peil
%peil_prim_j  % jaar gemiddelde peil
%peil_prim_w  % winterpeil
%peil_prim_z  % zomerpeil
%oppervlakte_owp
cdr_prim=cdr_prim_default*ones(size(peil_prim_j));  % default by lack of NHI data
IRIV=find(~isnan(bodemhoogte_prim_j) & ~isnan(peil_prim_j));
RIVindices=cellIndices(IRIV,[NROW,NCOL,NLAY]);
%plot(xm(RIVindices(:,2)),ym(RIVindices(:,1)),'.');
RIV=[RIVindices(:,[3,1,2]),peil_prim_j(IRIV),...
    oppervlakte_owp(IRIV)./cdr_prim(IRIV),bodemhoogte_prim_j(IRIV)];

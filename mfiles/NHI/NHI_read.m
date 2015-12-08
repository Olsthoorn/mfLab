function NHI_read(ziplist,URL,xLim,yLim)
%NHI_READ  reads NHI datafiles and save in mat file to minimize space
%
% Example:
%   NHI_READ(pth [,xLim,yLim])  looks for ASC files in the path and
%
%   reads them in one by one, selecting the part accordign to xLim, yLim
%   and saves to NHI_IX1_IX2_IY1_IY2.mat or to nhi struct if spcicified
%
% ToDo: needs thorough check and consistent link with other NHI files
% (130429)
%
% JB 090722 TO 110425

fprintf('reading NHI ASC data files dir ''%s''\n',home);
fprintf('%s\n',datestr(now));

%% Initialize size of AGV model
if nargin<4
    xLim=[100000 155000];
    yLim=[450000 500000];
end

if nargin<3 | ~exist('URL','var') | isempty(URL)
    URL='http://www.NHI.nu/Downloads';
end

for i=1:length(ziplist)


% =================== READ ASC FILES ===================== 
fprintf('Reading ASC files from path ''%s''\n',pth);

d=dir([pth '*.asc']);  % ascii NHI data set files

%% Determine size model -- only red the file headers

NLAY=0;
for i=1:length(d)
    [ptype,layNr]=identify(d(i).name);
    if ptype=='k' && layNr>NLAY,
        NLAY=layNr;
        if layNr==1
            [dummy,xGr,yGr,xm,ym,Ix,Iy]=NHI_readASC(d(i).name,pth,xLim,yLim);
            NCOL=length(xm);
            NROW=length(ym);
        end
    end
end  

fprintf('AGV model: NCOL,NROW,NLAY= %d %d %d\n',NCOL,NROW,NLAY);

%% Allocate space for matrices
KD    =NaN(NROW,NCOL,NLAY);
VCONT =NaN(NROW,NCOL,NLAY-1);
THK   =NaN(NROW,NCOL,NLAY);
IBOUND=NaN(NROW,NCOL,NLAY);
STRTHD=NaN(NROW,NCOL,NLAY);

%% Read and fill matrices

for i=1:length(d)
    FName=d(i).name;
    [ptype,layNr]=identify(FName);
    switch ptype
        case 'k'
            fprintf('reading kD layer %d from file ''%s''\n',layNr,FName);
            KD(:,:,layNr)=readASC(FName,pth,xLim,yLim);
        case 'v'
            fprintf('reading VCONT layer %d from file ''%s''\n',layNr,FName);
            VCONT(:,:,layNr)=readASC(FName,pth,xLim,yLim);
        case 'n'
            fprintf('reading STRTHD layer %d from file ''%s''\n',layNr,FName);
            STRTHD(:,:,layNr)=readASC(FName,pth,xLim,yLim);
        case 'd'
            fprintf('reading THK layer %d from file ''%s''\n',layNr,FName);
            THK(:,:,layNr)=readASC(FName,pth,xLim,yLim);
        case 'i'
             fprintf('reading IBOUND layer %d from file ''%s''\n',layNr,FName);
             IBOUND(:,:,layNr)=readASC(FName,pth,xLim,yLim);
        case 'f'
            fprintf('reading ANIFCT (anisotropy factor) from file ''%s''\n',FName);
            ANIFCT =readASC(FName,pth,xLim,yLim);
        case 'h'
            fprintf('reading HK (anistropy angle) from file ''%s''\n',FName);
            HK =readASC(FName,pth,xLim,yLim);
        case 'r'
            fprintf('reading RCH from file ''%s''\n',FName);
            RCH =readASC(FName,pth,xLim,yLim);
        case 'c'
            fprintf('reading CNC from file ''%s''\n',FName);
            %CNC=readASC(FName,pth,xLim,yLim);
        case 'C'
            fprintf('reading CONC layer %d from file ''%s''\n',layNr);
            %CONC(:,:,layNr)=readASC(FName,pth,xLim,yLim);
        case 'a'
            fprintf('reading AHN from file ''%s''\n',FName);
           AHN=readASC(FName,pth,xLim,yLim);
        otherwise
            error('Don''t know wat to do with this ptype %s',ptype);
    end
end

% JOS zoekt verder uit hoe we met IBOUND omgaan in het AGV model
Z=NaN(NROW,NCOL,NLAY+1);

Z(:,:,1)=AHN;

thkmin=0.1; % minimum thickness of layers


for iLay=1:NLAY
    Z(:,:,iLay+1)=Z(:,:,iLay)-max(thkmin,THK(:,:,iLay));
    if iLay>2, IBOUND(:,:,iLay)=IBOUND(:,:,2); end
end

IBOUND([1,end],:,:)=-1;
IBOUND(:,[1,end],:)=-1;
IBOUND(STRTHD<Z(:,:,2:end))=1;  % perhaps set on 0 (inactive)

%showprov([shppath shpfile]);
%contour(xm,ym,A)

% ================== READ SCD FILES =============================
pth=[base 'scd\'];
fprintf('Reading SCD files from path ''%s''\n',pth);

d=dir([pth '*.scd']);  % scd   NHI data set files

for i=1:length(d)
    FName=d(i).name;
    [ptype,layNr]=identify(FName);
    switch ptype
        case 'd'
            fprintf('reading SCD file ''%s'' for DRN\n',FName);
            DRN=readSCD(FName, pth,'d',Ix,Iy);
        case 'g'
            fprintf('reading SCD file ''%s'' for GHB\n',FName);
            GHB=readSCD(FName, pth,'g',Ix,Iy);
        case 'w'
           fprintf('reading SCD file ''%s'' for WEL\n',FName);
            WEL=readSCD(FName, pth,'w',Ix,Iy);
        case 'r'
           fprintf('reading SCD file ''%s'' for RIV\n',FName);
            RIV=readSCD(FName, pth,'r',Ix,Iy);
        otherwise
            fprintf('Skipping file ''%s'' due to ptype ''%s''\n',FName,ptype);
    end
end
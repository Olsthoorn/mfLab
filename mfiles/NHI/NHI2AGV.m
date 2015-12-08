%NHI2AGV extracts AGV model from NHI datafiles
%
% JB 090722

%base='D:\data\NHI\model\input\mf\v067\';
base='Z:\tolsthoorn On My Mac\GRWMODELS\MYWORK\AGVmodel\NHI2AGV\';
home=[base 'mfiles'];

path(path,home);

cd(home);
fprintf('Running NHI2AGV from dir ''%s''\n',home);
fprintf('%s\n',datestr(now));

%% Initialize size of AGV model
xLimAGV=[100000 155000];
yLimAGV=[450000 500000];

% =================== READ ASC FILES ===================== 
pth=[base 'asc\'];
fprintf('Reading ASC files from path ''%s''\n',pth);

d=dir([pth '*.asc']);  % ascii NHI data set files
%% Determine size of AGV model

NLAY=0;
for i=1:length(d)
    [ptype,layNr]=identify(d(i).name);
    if ptype=='k' && layNr>NLAY,
        NLAY=layNr;
        if layNr==1
            [dummy,xGr,yGr,xm,ym,IxLim,IyLim]=readASC(d(i).name,pth,xLimAGV,yLimAGV);
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
            KD(:,:,layNr)=readASC(FName,pth,xLimAGV,yLimAGV);
        case 'v'
            fprintf('reading VCONT layer %d from file ''%s''\n',layNr,FName);
            VCONT(:,:,layNr)=readASC(FName,pth,xLimAGV,yLimAGV);
        case 'n'
            fprintf('reading STRTHD layer %d from file ''%s''\n',layNr,FName);
            STRTHD(:,:,layNr)=readASC(FName,pth,xLimAGV,yLimAGV);
        case 'd'
            fprintf('reading THK layer %d from file ''%s''\n',layNr,FName);
            THK(:,:,layNr)=readASC(FName,pth,xLimAGV,yLimAGV);
        case 'i'
             fprintf('reading IBOUND layer %d from file ''%s''\n',layNr,FName);
             IBOUND(:,:,layNr)=readASC(FName,pth,xLimAGV,yLimAGV);
        case 'f'
            fprintf('reading ANIFCT (anisotropy factor) from file ''%s''\n',FName);
            ANIFCT =readASC(FName,pth,xLimAGV,yLimAGV);
        case 'h'
            fprintf('reading HK (anistropy angle) from file ''%s''\n',FName);
            HK =readASC(FName,pth,xLimAGV,yLimAGV);
        case 'r'
            fprintf('reading RCH from file ''%s''\n',FName);
            RCH =readASC(FName,pth,xLimAGV,yLimAGV);
        case 'c'
            fprintf('reading CNC from file ''%s''\n',FName);
            %CNC=readASC(FName,pth,xLimAGV,yLimAGV);
        case 'C'
            fprintf('reading CONC layer %d from file ''%s''\n',layNr);
            %CONC(:,:,layNr)=readASC(FName,pth,xLimAGV,yLimAGV);
        case 'a'
            fprintf('reading AHN from file ''%s''\n',FName);
           AHN=readASC(FName,pth,xLimAGV,yLimAGV);
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
            DRN=readSCD(FName, pth,'d',IxLim,IyLim);
        case 'g'
            fprintf('reading SCD file ''%s'' for GHB\n',FName);
            GHB=readSCD(FName, pth,'g',IxLim,IyLim);
        case 'w'
           fprintf('reading SCD file ''%s'' for WEL\n',FName);
            WEL=readSCD(FName, pth,'w',IxLim,IyLim);
        case 'r'
           fprintf('reading SCD file ''%s'' for RIV\n',FName);
            RIV=readSCD(FName, pth,'r',IxLim,IyLim);
        otherwise
            fprintf('Skipping file ''%s'' due to ptype ''%s''\n',FName,ptype);
    end
end
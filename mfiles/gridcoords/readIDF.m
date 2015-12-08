function [out,incol,inrow]=readIDF(fname,option,arg3,fmt,fac)
%READIDF reads IMOD's IDF file (Peter Vermeulen, Deltares)
%
% Example:
%    [out,incol,inrow] = readIDF(fname,option,arg3,fmt,fac) -- reads MIPWA idf-file
% 
% fname:      fname with entire MIPWA grid
% option:     'cs'
% arg3:       [x(:) y(:)] cross section
%
% option:     'window'
% arg3:       [xmin xmax ymin ymax] curout
%
% fmt         format to convert contents of file to (float); default 'int32'
% fac         factor by which values are to be multipllied (default 100)
%
% Copyright 2010 Frans Schaars, Artesia

if nargin==0
    % fname='..\data\mipwa\MODELLAGEN\top_laag1.IDF';
    % fname='..\data\mipwa\kd_new\KD_NEW1.IDF';
    %     option='cs';
    %     xy =[ 252447      553001
    %           252500      553036
    %           252508      553101];
    option='window';
    axlim=[   210000      290000      530000      620000];
    xy =[ 251447      580132
        244342      576447
        241711      572237
        239079      566184
        234605      559868
        232500      553026
        230658      549079
        227763      545921
        225921      541447
        222237      536184];

elseif nargin==1
    option='all';
elseif strcmp(option,'window')
    axlim=arg3;
elseif strcmp(option,'cs')
    xy=arg3;
elseif strcmp(option,'csint')
    xy=arg3;
end
if nargin<4
    fmt='int32';
    fac=100;
end



fid  = fopen(fname);
idum = fread(fid,1,'int');
ncol = fread(fid,1,'int');
nrow = fread(fid,1,'int');
xmin = fread(fid,1,'float');
xmax = fread(fid,1,'float');
ymin = fread(fid,1,'float');
ymax = fread(fid,1,'float');
amin = fread(fid,1,'float');
amax = fread(fid,1,'float');
xnodata = fread(fid,1,'float');

regular = fread(fid,1,'int');
if regular==0
    dxy=fread(fid,1,'float');
    x=xmin+0.5*dxy:dxy:xmax-0.5*dxy;
    y=ymin+0.5*dxy:dxy:ymax-0.5*dxy;
else
    error('regular is 1');
end
% randen
% xr=[x-0.5*dxy x(end)+0.5*dxy];
% yr=[y-0.5*dxy y(end)+0.5*dxy];


switch option
    case 'all'
        f=ones(nrow,ncol,fmt);
        nrowsub=1;
        tmp=NaN(nrowsub,ncol);
        for irow=1:nrow/nrowsub
            erin=1+nrowsub*(irow-1):nrowsub*irow;
            erin(erin>nrow)=[];
            if ~isempty(erin)
                tmp(:)=fread(fid,[ncol,length(erin)],'float')';
                f(nrow+1-erin,:)=(tmp*fac);
            end
        end
        out.x=x;
        out.y=y;
        out.v=f;
    case 'window'
        xmin=axlim(1);
        xmax=axlim(2);
        ymin=axlim(3);
        ymax=axlim(4);
        inrow=find(y>=ymin&y<=ymax);
        incol=find(x>=xmin&x<=xmax);
        f=ones(length(inrow),length(incol),fmt);
        skip=4*(nrow-max(inrow))*ncol;
        fseek(fid,skip,'cof');
        for irow=1:length(inrow)
            tmp=fread(fid,ncol,'float')';
            f(length(inrow)+1-irow,:)=(tmp(incol)*fac);
        end
        x=x(incol);
        y=y(inrow);
        out.x=x;
        out.y=y;
        out.v=f;
    case 'csint'

        % bepaal de waarden van de opgegeven doorsnede
        xp=xy(:,1);
        yp=xy(:,2);
        % hoekpunten van een grid dat groot genoeg is
        xmin=min(xp)-2*dxy;
        xmax=max(xp)+2*dxy;
        ymin=min(yp)-2*dxy;
        ymax=max(yp)+2*dxy;

        % zoek een window waarin de doorsnede past.
        inrow=find(y>=ymin&y<=ymax);
        incol=find(x>=xmin&x<=xmax);

        % haal de waarden voor dit window uit de file
        f=ones(length(inrow),length(incol),fmt);
        skip=4*(nrow-max(inrow))*ncol;
        fseek(fid,skip,'cof');
        for irow=1:length(inrow)
            tmp=fread(fid,ncol,'float')';
            f(length(inrow)+1-irow,:)=(tmp(incol)*fac);
        end

        % pas x en y (middens van de cellen) aan voor dit window
        x=x(incol);
        y=y(inrow);

        [x_int,y_int]=interpLine(xp,yp,1,'dist',1);
        v=interp2(x,y,single(f),x_int,y_int,'nearest');
        out.x=x_int;
        out.y=y_int;
        out.v=v;

    case 'cs'

        % bepaal de waarden van de opgegeven doorsnede
        xp=xy(:,1);
        yp=xy(:,2);
        % hoekpunten van een grid dat groot genoeg is
        xmin=min(xp)-2*dxy;
        xmax=max(xp)+2*dxy;
        ymin=min(yp)-2*dxy;
        ymax=max(yp)+2*dxy;

        % zoek een window waarin de doorsnede past.
        inrow=find(y>=ymin&y<=ymax);
        incol=find(x>=xmin&x<=xmax);

        % haal de waarden voor dit window uit de file
        f=ones(length(inrow),length(incol),fmt);
        skip=4*(nrow-max(inrow))*ncol;
        fseek(fid,skip,'cof');
        for irow=1:length(inrow)
            tmp=fread(fid,ncol,'float')';
            f(length(inrow)+1-irow,:)=(tmp(incol)*fac);
        end

        % pas x en y (middens van de cellen) aan voor dit window
        x=x(incol);
        y=y(inrow);
        out=cs(xp,yp,x,y,f);
end
fclose(fid);

%%
if nargout==0
    figure
    if strcmp(option,'all')
        imagesc(x(round(end/2):end),y(round(end/2):end),f(round(end/2):end,round(end/2):end))
    else
        imagesc(x,y,f);
    end
    axis image
    axis xy
    caxis([-32768       16300]);colorbar


    axlim=[   210000      290000      530000      620000];
    axis(axlim)
    %%
    figure
    %[X,Y]=meshgrid(x+.5*dxy,y+.5*dxy);
    %mesh(X,Y,0*Y)
    imagesc(x,y,f)
    colorbar
    axis image
    view(2)
    hold on
    
%     plot(xx,yy,'marker','.')
%     axis auto
%     plot(xp,yp,'o')
%     %    text(xx,yy,num2str((1:length(xx))'))
%     text(xx,yy,num2str(v(:)))

end

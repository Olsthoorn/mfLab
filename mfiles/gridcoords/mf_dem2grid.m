function newdem=mf_dem2grid(dem,eGr,nGr)
%MF_DEM2GRID make a dem using grid coordinates
%
% Example:
%    newdem=mf_dem2grid(dem,eGr,nGr)
%
%    newdem has same structure as dem
%    dem is a struct with fields as obtained from
%    mf_getdemfromtiff(TiffFileName)
%
% See also mf_getdemfromtiff
%
% ToDo: generalize and turn dem into a useful object
%
% TO 110530

z=dem.z; dem.z=[];

newdem.eGr=eGr; % easting (longitude)
newdem.nGr=nGr; % northing (latitude)

newdem.em=0.5*(newdem.eGr(1:end-1)+newdem.eGr(2:end)); % longitude
newdem.nm=0.5*(newdem.nGr(1:end-1)+newdem.nGr(2:end)); % latitude

newdem.BoundingBox(:,1)=newdem.eGr([1 end]);
newdem.BoundingBox(:,2)=newdem.nGr([1 end]);

newdem.Width =length(newdem.eGr)-1;
newdem.Height=length(newdem.nGr)-1;

newdem.z=NaN(newdem.Height,newdem.Width);

tic
Ie1=NaN(1,newdem.Width);
Ie2=NaN(1,newdem.Width);
EGr=dem.eGr;
for ie=1:newdem.Width
    Ie1(ie)=find(EGr>=eGr(ie  ),1,'first');
    Ie2(ie)=find(EGr<=eGr(ie+1),1,'last');
end
toc
tic

In1=NaN(1,newdem.Height);
In2=NaN(1,newdem.Height);
NGr=dem.nGr;
for in=1:newdem.Height
    In1(in)=find(NGr>=nGr(in  ),1,'first');
    In2(in)=find(NGr<=nGr(in+1),1,'last');
end
toc

tic
for ie=1:newdem.Width
    zz=mean(z(:,Ie1(ie):Ie2(ie)),2);
    zzz=NaN(1,newdem.Height);
    for in=1:newdem.Height
        zzz(in)=mean(zz(In1(in):In2(in)));
    end
    newdem.z(:,ie)=zzz';
end
toc

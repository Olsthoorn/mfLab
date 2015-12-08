function newdem=mf_demCoarse(dem,nE,nN)
%MF_DEMCOARSE generates a coarser dem (digital elevation model)
%
% Example:
%    newdem=mf_demCoarse(olddem,nE,nN)
%
%    newdem has same structure as dem but is nE,nN coarser
%    dem is a struct with fields as obtained from
%
% See also: geotiffinfo geotiffread
%
% ToDo: generalize and turn into a useful object
%
% TO 110530

Z=dem.z; dem.z=[];

newdem.PixelScale=dem.PixelScale.*[nE nN 1]';

newdem.eGr=dem.BoundingBox(1,1):newdem.PixelScale(1):dem.BoundingBox(2,1);
newdem.nGr=dem.BoundingBox(1,2):newdem.PixelScale(2):dem.BoundingBox(2,2);

newdem.em=0.5*(newdem.eGr(1:end-1)+newdem.eGr(2:end));
newdem.nm=0.5*(newdem.nGr(1:end-1)+newdem.nGr(2:end));
newdem.BoundingBox(:,1)=newdem.eGr([1 end]);
newdem.BoundingBox(:,2)=newdem.nGr([1 end]);

newdem.Width =length(newdem.eGr)-1;
newdem.Height=length(newdem.nGr)-1;

newdem.z=ones(newdem.Width,newdem.Height);

tic
for ie=1:newdem.Width-1
    IX=(ie-1)*nE+(1:nE);
    zz=mean(Z(:,IX),2);
    for in=1:newdem.Height-1
        IY=(in-1)*nN+(1:nN);
        newdem.z(in,ie)=mean(zz(IY));
    end
end
toc

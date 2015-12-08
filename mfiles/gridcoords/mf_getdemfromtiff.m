function dem=mf_getdemfromtiff(FName)
%MF_GETDEMFROMTIFF read a DEM from a TIFF file
%
% Example:
%   dem=mf_getdemfromtiff(TiffFilenName)
%
% ToDo: check and generalize by converting to a useful object (TO 130428)
%
% See also: mf_cleandem, mf_dem2grid mf_demCoarse
%
% TO 110530

dem=geotiffinfo(FName);
dem.z=geotiffread(FName);
dem.z=flipud(dem.z);

dem.eGr=dem.BoundingBox(1,1):dem.PixelScale(1):(dem.BoundingBox(2,1)+dem.PixelScale(1)/100);
dem.nGr=dem.BoundingBox(1,2):dem.PixelScale(2):(dem.BoundingBox(2,2)+dem.PixelScale(1)/100);

dem.em=0.5*(dem.eGr(1:end-1)+dem.eGr(2:end));
dem.nm=0.5*(dem.nGr(1:end-1)+dem.nGr(2:end));

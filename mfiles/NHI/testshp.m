%TESTSHP tests shape file
%
% Example:
%    testshp
%
% See also: readshp
%
% TO 090728

FName='Z:\tolsthoorn On My Mac\GRWMODELS\NHI_and_AGV\Bibliotheek_NHI\districten\districten';

[shape,SHAPE,data]=readshp(FName);
plotshp(shape)
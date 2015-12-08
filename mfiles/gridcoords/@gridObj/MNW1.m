function well = MNW1(gr,basename,HK,wellOrWellSheetNm,varargin)
% well = gr.well(basename,HK,wellOrWellSheetName [,QsheetNm [,CsheetNm [,use]
%
%  calls 
%  well = MNW1Obj(basename,well,gr,HK,QsheetNm,CsheetNm,use);
%
%  see wellObj
%
% SEE ALSO: gridObj/setWell wellObj/WEl wellObj/PNTSRC mfSetWells
%
%   TO 110426 120103 120408 121118
%
% Copyright 2009-2012 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later


well = MNW1Obj(basename,wellOrWellSheetNm,gr,HK,varargin{:});




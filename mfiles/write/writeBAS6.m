function writeBAS6(basename,bas)
%WRITEBAS6 writes input file for MODFLOW's basic (BAS6) package
%
% Example:
%    writeBAS6(basename,bas)
%
% TO 070630


% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen([basename,'.',bas.ext],'wt');

%0.
fprintf(fid,'# MATLAB  writeBAS6 %s\n',datestr(now));
fprintf(    '# MODFLOW writeBAS6 %s\n',datestr(now));

%1.  optional words XSECTION CHTOCH FREE PRINTTIME SHOWPROGRESS
OPTIONS={};
if bas.FREE,   OPTIONS{1}='FREE'; end
if bas.CHTOCH, OPTIONS{2}='CHTOCH'; end

if ~isempty(OPTIONS)
     fprintf(fid,' %s',OPTIONS{:});
end
fprintf(fid,'\n');

ctrlRec = true;

%2.  if <0 constant head if > 0 compute head if 0 inactive (or lake in lake package)
for ilay=1:bas.NLAY
    warray(fid,bas.IBOUND(:,:,ilay),bas.unit,'(25I3)',sprintf('IBOUND{%d}',ilay),ctrlRec,bas.FREE);
end

%3. value used for inactive cells
fprintf(fid,'%10g     HNOFLO\n',bas.HNOFLO);

%4. initial heads
if numel(bas.STRTHD)==1
    for ilay=1:bas.NLAY
        warray(fid,bas.STRTHD          ,bas.unit,'(10E13.5)',sprintf('STRTHD{%d}',ilay),ctrlRec,bas.FREE);
    end
else
    for ilay=1:bas.NLAY
        warray(fid,bas.STRTHD(:,:,ilay),bas.unit,'(10E13.5)',sprintf('STRTHD{%d}',ilay),ctrlRec,bas.FREE);
    end
end
fclose(fid);

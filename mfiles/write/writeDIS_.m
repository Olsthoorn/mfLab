function writeDIS(basename,dis)
%WRITEDIS writes input file for MODFLOW's dicretization package (DIS)
%
% Example:
%    writeDIS(basename,dis)  --- writing discretization file
%
% TO 070630 100827 120426

% Copyright 2009-2012 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fid=fopen([basename,'.',dis.ext],'wt');

%0.
fprintf(fid,'# MATLAB  writeDIS %s\n',datestr(now));
fprintf(    '# MODFLOW writeDIS %s\n',datestr(now));

%1.
fprintf(fid,'%10d%10d%10d%10d%10d%10d      %s\n',...
    dis.GRID.Nlay,dis.GRID.Ny,dis.GRID.Nx,dis.NPER,dis.ITMUNI,dis.LENUNI,...
    'NLAY, NROW NCOL NPER ITMUNI LENUNI');

%2  Resistance at bottom of layer?? -- Quasi 3D flag for each layer
warray(fid,dis.GRID.LAYCBD(:)',dis.unit,'(40I2)','LAYCBD',false);

%3.
warray(fid,dis.GRID.dx(:)',dis.unit,'(10E12.3)','DELX or DELR',true,dis.FREE);

%4.
warray(fid,dis.GRID.dy(:)',dis.unit,'(10E12.3)','DELY or DELC',true,dis.FREE);

%5
warray(fid,dis.GRID.Z(:,:,1),dis.unit,'(10E15.6)','TOP of model',true,dis.FREE);

%6.
for i=2:dis.GRID.Nz+1
    warray(fid,dis.GRID.Z(:,:,i),dis.unit,'(10E15.6)',sprintf('BOT{%d}',i),true,dis.FREE);
end

%% FOR EACH STRESS PERIOD
%7. 
for iP=1:dis.NPER
%    fprintf(fid,'%10s%10d%10g',numsqueeze(dis.PERLEN(iP),'%g',10),dis.NSTP(iP),dis.TSMULT(iP));
    fprintf(fid,'%12.6g %10d %10g',dis.PERLEN(iP),dis.NSTP(iP),dis.TSMULT(iP));
    if dis.isTran(iP),
        fprintf(fid,'     TR');
    else
        fprintf(fid,'     SS');
    end;
    fprintf(fid,'     PERLEN NSTP TSMULT SS/TR\n');
end

fclose(fid);

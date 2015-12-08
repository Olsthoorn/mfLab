function  endp=readENDP(fname)
%READENDP reads end points file produced by MODPATH6
%
% Example:
%    endp=readENDP(fname);
%
% TO 120330

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%0.
fprintf('# MATLAB readENDP %s\n',datestr(now));

fid=fopen(fname,'r');

% CHARACTER*80 TEXT
endp.hdr=scanf(fid,'%80s',1);

% TREF
endp.tref=fread(fp,1,'float32');

p1=ftell(fp); contrRec(fp,bytes); reclen=ftell(fp)-p1;
fseek(fp,0,1); NP=(ftell(fp)-p1)/reclen;

if floor(NP)~=NP
    error('mfLab:readENDP:reclen',...
        'ReadENDP: Noninteger number of records(%f) of lengt (%d) in file %s\n',NP,reclen,fname)
else
    fprintf('%d endpoint records found in binary modpath endpoint file.\n',NP)
end

[ep(NP).iz_e,...
    ep(NP).idx_e,...
    ep(NP).ex,...
    ep(NP).ey,...
    ep(NP).ez,...
    ep(NP).et,...
    ep(NP).bx,...
    ep(NP).by,...
    ep(NP).bz,...
    ep(NP).idx_b,...
    ep(NP).iz_b,...
    ep(NP).itm,...
    ep(NP).ipc,...
    ep(NP).t_release] = contrRec(fp,bytes);

for i=1:NP
[ep(i).iz_e,...
    ep(i).idx_e,...
    ep(i).ex,...
    ep(i).ey,...
    ep(i).ez,...
    ep(i).et,...
    ep(i).bx,...
    ep(i).by,...
    ep(i).bz,...
    ep(i).idx_b,...
    ep(i).iz_b,...
    ep(i).itm,...
    ep(i).ipc,...
    ep(i).t_release] = contrRec(fp,bytes);
end
    

function [iz_e,idx_e,ex,ey,ez,et,bx,by,bz,idx_b,iz_b,itm,ipc,t_release]=contrRec(fp)
% [iz,idx,ex,ey,ez,et,bx,by,bz,idx,izb,itm,ipc]=contrRec(fp,bytes)
% --- reads a complete layer from a MODPATH end points file binary file
% TO 12-330

% Particle data records contain the following items:
% 1. Zone code for the cell containing the final location of the particle
iz_e  = fread(fp, 1,'int32'  );

% 2. Global node number for the cell containing the final location
idx_e = fread(fp, 1,'int32');

% 3. Global coordinate in the x-direction (J index direction) for the final location
ex  = fread(fp, 1,'float32');

% 4. Global coordinate in the y-direction (I index direction) for the final location
ey  = fread(fp, 1,'float32');

% 5. Local coordinate for the z-direction within the grid cell
ez  = fread(fp, 1,'float32');

% 6. Total tracking time
et  = fread(fp, 1,'float32');

% 7. Global coordinate in the x-direction for starting location
bx  = fread(fp, 1,'float32');

% 8. Global coordinate in the y-direction for starting location
by  = fread(fp, 1,'float32');

% 9. Local coordinate in the z-direction within the cell for starting location
bz  = fread(fp, 1,'float32');

% 10. Global node number for the cell containing starting location
idx_b = fread(fp, 1,'int32'  );

% 11. Zone code for cell containing starting location
iz_b = fread(fp, 1,'int32'  );

% 12. Cumulative MODFLOW time step number corresponding to the time of release
itm = fread(fp, 1,'int32'  );

% 13. Particle termination code, IPCODE
ipc = fread(fp, 1,'int32'  );

%14. Release time
t_release = fread(fp, 1,'float32');

fclose(fid);

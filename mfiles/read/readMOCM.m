%function mocm = readMOCM(fname,pth)
%READMOCM reads the mocmain.dat file which is used to run MocDens3D
%
% Example:
%
% ToDo: currently a hard wired function. Generalize and test.
%
% TO 090806

%clear all
%close all

pth='C:\Users\David\Desktop\GRWMODELS\MYWORK\mocdens3d\';
fname='mocmain.dat';

fprintf('# MATLAB readMOC %s\n',datestr(now));

fid=fopen([pth fname],'r');
skipmodflowcomments(fid);

mocm.name = fgets(fid);

mocm.size = fscanf(fid,'%f',23);
mocm.NLAY = mocm.size(2);
mocm.NROW = mocm.size(4);
mocm.NCOL = mocm.size(6);

mocm.DENS = NaN(mocm.NROW,mocm.NCOL,mocm.NLAY);
fprintf('Reading the density values.\n');
for iLay = 1:mocm.NLAY
    mocm.DENS(:,:,iLay) = mudread(fid,[mocm.NROW,mocm.NCOL]);
end

mocm.afterDENS = fscanf(fid,'%f',10);

mocm.SOURCES = NaN(mocm.NROW,mocm.NCOL,mocm.NLAY);
fprintf('Reading the sources data values.\n');
for iLay = 1:mocm.NLAY
    mocm.SOURCES(:,:,iLay) = mudread(fid,[mocm.NROW,mocm.NCOL]);
end

s = fgets(fid);
fprintf('Reading the Longitudinal Dispersivity.\n');
mocm.LD = fscanf(fid,'%f',mocm.NLAY);
s = fgets(fid); s = fgets(fid);
fprintf('Reading the Horizontal Transverse Dispersivity.\n');
mocm.HTD = fscanf(fid,'%f',mocm.NLAY);
s = fgets(fid); s = fgets(fid);
fprintf('Reading the Vertical Transverse Dispersivity.\n');
mocm.VTD = fscanf(fid,'%f',mocm.NLAY);
s = fgets(fid); s = fgets(fid);
fprintf('Reading the Retardation Factor.\n');
mocm.RF = fscanf(fid,'%f',mocm.NLAY);


s = fgets(fid);
mocm.THICKNESS = NaN(mocm.NROW,mocm.NCOL,mocm.NLAY);
mocm.PRSTY = NaN(mocm.NROW,mocm.NCOL,mocm.NLAY);
fprintf('Reading the layers thickness values.\n');
fprintf('Reading the layers porosity values.\n');
for iLay = 1:mocm.NLAY
    mocm.THICKNESS(:,:,iLay) = mudread(fid,[mocm.NROW,mocm.NCOL]);
    mocm.PRSTY(:,:,iLay) = mudread(fid,[mocm.NROW,mocm.NCOL]);
end

save('mocm.mat','mocm','BOTTOM','TOP','NLAY','NSURF','Nx','Ny');

fclose(fid);

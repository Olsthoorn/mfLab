function readADF
%READADF potentially reads all grid datasets of the directory one by one.
%
% transforms the data to an ESRI ASCII grid which is further transformed by
% Matlab to a mat file of the same name.
%
% A complete transfer of all files is an immense task which  will take
% hours to days perhaps

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fprintf('#readADF  %s\n',datestr(now));


P='Z:\tolsthoorn On My Mac\GRWMODELS\TNO-REGIS\mirone131\';
G='Z:\tolsthoorn On My Mac\GRWMODELS\TNO-REGIS\amstel gooi en vechtstreek\hydrogeologie\grids\';

fprintf('Reading from directory\n%s\n',P);

cd(G)

d=dir([G '*']); d([1 2])=[];

fprintf('%d grid datasets\n',sum([d.isdir]));

fplog=fopen('readADF.log','a');

for i=1:length(d)
    if d(i).isdir && ~exist([d(i).name,'.mat'],'file')
         system(['"', P, 'gdal_translate.exe" -of AAIGrid "',d(i).name,...
             '" "',d(i).name,'-ASCII"']);
         try
             A=readASC([d(i).name,'-ASCII']);
             save([d(i).name,'.mat'],'A');
             delete([d(i).name,'-ASCII']);
             fprintf(fplog,'%3d: %-20s: ok:    %s\n',i,d(i).name,datestr(now));
             fprintf(      '%3d: %-20s: ok:    %s\n',i,d(i).name,datestr(now));
         catch ME
             fprintf(fplog,'%3d: %-20s: error: %s\n',i,d(i).name,ME.message);
             fprintf(      '%3d: %-20s: error: %s\n',i,d(i).name,ME.message);
         end
    else
        fprintf(fplog,'%3d: %-20s: skipped\n',i,d(i).name);
        fprintf(      '%3d: %-20s: skipped\n',i,d(i).name);
    end
end

fclose('all');
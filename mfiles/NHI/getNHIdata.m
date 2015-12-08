%GETNHIDATA script to read the input files of the Netherlands Hydrologic Instrument (NHI)
%
% Example:
%   getNHIdata
%
% The files are read from
%    www.nhi.nu/downloads
%
% The names of the zip files are in workbook NHI.xls. They have been
% obtained by clicking op de links on the mentioned site and copying this
% link to the workbook. So this link may have to be updates when developers
% of the NHI change their minds.
%
% TO 120427

basename = 'NHI';
site     = ['http://www.' basename '.nu/'];

home   =['/Users/Theo/GRWMODELS/mflab/examples/CT5440/' basename '/'];
dlzips =[home 'Downloads1/'];
dlfiles=[home 'Downloads2/'];

fprintf('The data will be read form this site:\n%s\n', ...
         [site 'downloads']);

% Get lcoations form workbook (put in there by hand)
[nams,vals,txthdr,txtvals]=getExcelData([home basename],'variables','H');

%% 
j=strmatchi('Location',txthdr);

% download the zips
for i=1:length(txtvals),
    fprintf('Downloading ''%s'' ...',txtvals{i,j});
    [P F E]=fileparts(txtvals{i,j});
    try
        urlwrite(txtvals{i,j},[dlzips F E]);
        fprintf('done.\n');
    catch ME
        fprintf('%s\n',ME.message);
        fprintf('not done\n');
    end

end


%%
cd(dlfiles);
for i=1:length(txtvals)
    fprintf('unzipping ''%s'' ... ',txtvals{i,j});
    [P F E] = fileparts(txtvals{i,j});
    try
        unzip([dlzips F E],dlfiles);
        fprintf('done\n');
    catch ME
        fprintf('%s\n',ME.message);
        fprintf('not done\n');
    end
end

function D = mf_checkdir(varargin)
%MF_CHECKDIR verifies that output has been generated, i.e. the model run terminated normally.
%
% USAGE:
%    mf_checkdir()
%
% TO 120421

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

load name

fprintf(['Checking directory to see that model output files <<%s.*>>\n',...
          'are present and have a time stamp consistent with their input.\n'],basename);
     

d=[ dir([basename '.*']); dir('MT300*.UNC') ];

[~, I] = sort(vertcat(d.datenum),'descend');

d=d(I);

ext=cell(size(d));

for i=1:length(d)
    [~,~,ext{i}]=fileparts(d(i).name);
end

if nargin==0
    checklist={'.HDS','.BGT','.UCN'};
else
    checklist = varargin;
end

ihalf = round(length(d)/2);

err=[];

for i=1:length(checklist)
    fl = dir(['*' checklist{i}]);
    if isempty(fl)
        err{end+1}=sprintf('File <<*%s>> not found\n',checklist{i}); %#ok
    else
        j = strmatchi(fl(1).name,{d.name});
        if j(1)>ihalf
            err{end+1} = sprintf(['File <<%s>> not up to date (older than input)\n',...
                   'See that OC package in on in the NAM sheet\n',...
                   'See that IHDDFL IBUDFL ICBCFL in the PER sheet have the correct\n',...
                   'frequency and see that Hdpr Ddpr Hdsv Ddsv are on in the PER sheet,\n',...
                   'at least for some of them as desired.\n',...
                   'If all these parameters are set correctly, then check the LST file\n',...
                   'to see if the model has normal termination or chrashed somewhere along the way.'],fl.name); %#ok
        end
    end
end
        
if ~isempty(err)
    D=show(d);
    warning('mfLab:mf_checkdir:notUpToDataOrFileMissing',...
        [err{:}]);
else
    fprintf('...directory ok\n');
end

end

function D = show(d)
%D = show(d) shows the info on the basename and UCN files in this directory

    D=cell(length(d),4);

    fprintf('\n%25s%25s%12s%25s\n','Name', 'Date', 'Bytes', 'Date\n');

    for id=1:length(d)
            fprintf('%25s',d(id).name);
            fprintf('%25s',d(id).date);
            fprintf('%12d',d(id).bytes);
            fprintf('%25s',datestr(d(id).datenum));
            fprintf('\n');

            D(id,:) = {d(id).name, d(id).date, d(id).bytes,d(id).datenum};            
    end

end

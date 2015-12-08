function wel=mf_getwells(basename,welsheetnm)
% MF_GETWELLS: Collect well info from wellsheetnm into wel struct
%
% USAGE:
%    wel=mf_getwells(basename,welsheetnm)
%
%    wel is a wel sruct with all relevant information about the wel
%    welsheetnmae is the name of the sheet in the excel file basename
%
% SEE ALSO: mf_conf, mf_zone mf_plotConf
%
%   TO 111012

[welnams,welvals,weltxthdr,weltxt]=getExcelData(basename,welsheetnm,'hor');

for iw=size(welvals,1):-1:1 % running backward to allow removing while checking
    
    for j=1:length(weltxthdr)
        wel(iw).(weltxthdr{j})=weltxt{iw,j};
    end
    for j=1:length(welnams)
        wel(iw).(welnams{j})=welvals(iw,strmatchi(welnams{j},welnams,'exact'));
    end
end

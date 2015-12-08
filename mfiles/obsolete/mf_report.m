function mf_report(what,basename)
% mf_report(what,basename)   simple reporter about information used from worksheet
% what is one of {'MFLOW', 'MT3D' 'SEAWAT' 'SWI' 'NAM'}
% requires update and extension to all other worksheets in workbook.
% requires specific code per worksheet, that includes the package connected  to the parameter,
% which is currenty ignored
% TO 120415

fprintf('\n%s: Reporting on <<%s>> using excel file <<%s>>,mdl,basename:\n',mfilename,what,basename);

switch lower(what)
    case {'mflow', 'modflow' 'mf'}
        dir='v';
        [nams,vals]=getExcelData(basename,'MFLOW',dir);
    case {'mt3d', 'mt3dms', 'mt'}
        dir='v';
        [nams,vals]=getExcelData(basename,'MT3D',dir);
    case {'swt', 'seawat', 'swt_v4'}
        dir='v';
        [nams,vals]=getExcelData(basename,'SEAWAT',dir);
    case {'swi'}
        dir='v';
        [nams,vals]=getExcelData(basename,'MFLOW',dir);
    case {'nam'}
        dir='v';
        [nams,vals]=getExcelData(basename,'NAM',dir);
    otherwise
        error('mf_report,unknown option %s',what);
end

if dir=='v'
    for j=1:length(nams)
        fprintf('%15s',nams{j});
        for i=1:size(vals,2)
            fprintf('%12g',vals(j,i));
        end
        nl;
    end
else
    for iCol=1:length(nams), fprintf('%12s',nams{iCol}); end; nl;
    for iRow=1:size(vals,1)
        for iCol=1:length(nams),
            fprintf('%12g',vals(iCol));
        end
        nl;
    end
end
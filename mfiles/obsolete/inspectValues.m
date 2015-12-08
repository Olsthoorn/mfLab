function C = inspectValues(C,varargin)
% C = inspectValues(C,STCONC,field)
% inspect values in struct returned from readDat, readMT3D or readBud or a 3D array
% TO 130403

NCOL = 6;

[STCONC, varargin] = getNext(varargin,'double',[]);
[field , varargin] = getNext(varargin,'char','values');
[label , ~       ] = getNext(varargin,'char','');

if isfield(field,C)
    % ok
elseif isfield('term',C)
    field = 'term';
    if isempty(label)
        error(['%s: with field ''term'' you must specify a label to search for\n',...
                'For instance: ''FLOWRIGHTFACE'', ''FLOWLOWERFACE'', ''WELLS'', or any other\n',...
                'label in the budget-like file.'],mfilename);
    end
else
    % ok, try given field
end


if ~isempty(STCONC)
    C = deltaValues(C,STCONC);
end

statistics = NaN(numel(C),NCOL);

fprintf('Statistics of struct read by readData, readMT3D or readBud\n');
fprintf('Number of records = %04d\n',numel(C));

fprintf('Field name is: %s',field);
if ~isempty(label)
    fprintf(', label name is:  %s\n',label);
end
fprintf('\n');

fprintf('Size of field %s  = [%d %d %d]\n',size(C(1).(field)));
fprintf('\n');
fprintfs('%9s',{'it','max','mean','min','stdef','max-min'});
fprintf('\b');

if ~strcmp(field,'term')
    for i = 1:numel(C)
        statistics(i,:) = [
                i,             ...
        max(C(  i).(field)(:)) ...
        mean(C( i).(field)(:)) ...
        min(C(  i).(field)(:)) ...
        std(C(  i).(field)(:)) ...
        max(C(  i).(field)(:)) - min(C(i).(field)(:))];
    end
else
    for i = 1:numel(C)
        iL = strmatchi(label,C(i).label);
        if ~iL(1), continue; end

        iL=iL(1);
        statistics(i,:) = [ ...
            max(C(  i).(field){iL}(:)) ...
            mean(C( i).(field){iL}(:)) ...
            min(C(  i).(field){iL}(:)) ...
            stdev(C(i).(field){iL}(:)) ...
            max(C(  i).(field){iL}(:)) - min(C(i).(field){iL}(:))];
    end
end
display(statistics);

fprintfs('%9s',{'it','max','mean','min','stdef','max-min'});
fprintf('\n\n');

fprintfs('%11s',{'std(it)','std(max)','std(mean)','std(min)','std(std)','std(max-min)'});
fprintf('\n');
fprintf('%12g',std(statistics));
fprintf('\n');

msgId = 'mfLab:inspectValues:incompleteRecords';
if any(isnan(statistics(:,1)))
    warning('on',msgId);
    warning(msgId,'There are %d elements in the struct array that do not have your label <<%s>>.',...
        mfilename,sum(isnan(statistics(:,1))),label);
    warning('off',msgId);
end
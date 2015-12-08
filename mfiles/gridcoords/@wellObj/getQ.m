function [o,WEL] = getQ(o,basename,varargin)
%WELLOBJ/GETQ -- set the Q for wells
%
% USAGE
%  well = well.getQ(basename [,{sheetNm [prefix]},['index',indexFld])
%
%  EXAMPLE
%  well = well.getQ(basename);  % uses 'PER' as default sheetNm with 'Q'
%      prefix default column header in the PER sheet
%  well = well.getQ(basename,'PER');  % uses 'PER' as default sheetNm with 'Q'
%  well = well.getQ(basename,{'PER','Qday'});  % uses 'PER' as default sheetNm with 'Q'
%  well = well.getQ(basename,{'PER','Qday'},'index','Grp');  % uses 'PER' as default sheetNm with 'Q'
%     indexing of column to select based on field well.UserData.Grp instead of default well.nr.
%     this allows using other fields for selection of Q-columns.
%  well = well.getQ(basename,....,'fill',false)
%     prevents filling empty spreadsheet cells with values from above.
%
% basename = workbook name
% sheetNm = worksheet name in workbook containing, defaul 'PER'
% prefix is initial string of Q column header, default 'Q_'
%
% SEE ALSO: gridObj/setWell wellObj/WEl wellObj/PNTSRC mfSetWells
%
%   TO 110426 120103 120408 131119
%
% Copyright 2009-2012 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

[STATUS,SHEETS] = xlsfinfo(basename);
if isempty(STATUS)
    error('%s: Matlab can''t read excel workbook met basename <<%s>>',mfilename,basename);
end

[fill   ,varargin] = getProp(varargin,'fill',true);

% We may want to select a PER column for Q based on some UserData field,
% such as a group number. Here is how to specify it
[indexFld  ,varargin] = getProp(varargin,'index','');
index = ~isempty(indexFld);

[sheet_Q,varargin] = getNext(varargin,'char',[]);
[Sheet_Q,varargin] = getNext(varargin,'cell',{sheet_Q});

if ~isempty(Sheet_Q)
    I = strmatchi(Sheet_Q{1},SHEETS);
    if ~I(1)
        sheetNm = 'PER';
        prefix  = sheetNm;
    else
        sheetNm=Sheet_Q{1};
        Sheet_Q(1)=[];
        if ~isempty(Sheet_Q)
            prefix  = Sheet_Q{1};
        else
            prefix = 'Q';
        end
    end
else
    [Sheet_Q,varargin] = getNext(varargin,'char',[]);
    if isempty(Sheet_Q)
        error('%s: sheetNm must be of type char\n',mfilename);
    end
    I = strmatchi(Sheet_Q,SHEETS);
    if ~I(1)
        sheetNm = 'PER';
        prefix  = Sheet_Q;
    else
        sheetNm = Sheet_Q;
        prefix  = 'Q';
    end
end

if ~isempty(varargin)
    fprintf('%s: varargin list not empty after processing\n',mfilename);
end

[pernams,pervals,NPER]=getPeriods(basename,sheetNm,'fill',fill);

%% Add Data to the wells
if ~all(pervals(:,strmatchi('IHDDFL',pernams)) == pervals(:,strmatchi('ICBCFL',pernams)))
msgId = 'mfLab:wellObj_getQ:IHDDFLnotEqualICBCFL';
    warning(msgId,'on');    
    warning(msgId,...
        ['Not all output frequency flags IHDDFL equal ICBCFL in the PER sheet of workbook %s,\n',...
         'This implies that the output of heads|drawdown and that of the budget (Cell by Cell flows)\n',...
         'are not synchroneous. So be careful when matching these after reading\n',...
         'back in the binary heads|drawdowns and budget files.\n',...
         'You may use 9999 for both as a time step frequency flag to only get output after\n',...
         'every stress period and not after individual time steps.']);
    warning(msgId,'off');
end

Dt=pervals(:,strmatchi('PERLEN',pernams));

%% assert existance of prefix columns

QCOL = strmatchi(prefix,pernams);
if ~QCOL(1)
    error('%s: no flow columns with prefix <<%s>> exist in worksheet <<%s>> of workbook<<%s>>',...
        mfilename,prefix,sheetNm,basename);
end

% Attribute the correct flow column to each of the wells
for iw=numel(o):-1:1
    if index
        idx = o(iw).UserData.(indexFld);
    else
        idx = o(iw).nr;
    end

    qCol = strmatchi([prefix '_' num2str(idx)],pernams(QCOL),'exact');
    if ~qCol
        qCol= strmatchi([prefix  num2str(idx)] ,pernams(QCOL),'exact');
        if ~qCol
            qCol = strmatchi(prefix,pernams(QCOL),'exact');
            if ~qCol
                error(['%s: Prefix <<%s>> does not match any of the column headers in sheet <<%s>>\n',...
                       'neither is it unique. Therefore, your column request is incompatible with the headers.\n',...
                       'REMEDY: make sure that command\n',...
                       '        sprintf(%%s'',prefix%%d'',well(i).nr)\n',...
                       '        matches one of the headers in worksheet <<%s>> of workbook <<%s>>'],...
                        mfilename,prefix,sheetNm,sheetNm,basename);
            end
        end
    end
    qCol = QCOL(qCol); % number in the sheet

    o(iw).Q  = pervals(:,qCol)'; % as a horizontal vector for easy inspection
end

for iw=1:numel(o)
    % Get time step length
    o(iw).Dt = Dt';
    
    % Set o t. Note that this may be overruled by o.setCout(C,iComp).
    o(iw).t  = cumsum(o(iw).Dt,2);
    % Applications must verify wheather t starts at zero or at Dt !
end


if nargout<2, return; end


%% Nargout>1, so that the WEL array is requested
%  this array has [SP L R C Q] fields

% Count how many cells we have in total over all wells
% This requires that the wells have been set by o.toGrid(gr,HK)
% i.e. every o must have its grid information on board.

LRC = vertcat([o.LRC]);
NCell = size(LRC,1);

if isempty(LRC)
    error('mfLab:wellObj_getQ:noLRCdata',...
        '%s: the wells have not grid info on board, run o=o.toGrid(gr,HK) first.',mfile);
end

%% All wells carry their grid  info and fow info. So we can proceed with WEL

WEL = cell(NCell,1);
u   = ones(NPER,1);

k=0;
for iw=1:length(o)
    for i = 1:length(o(iw).idx)
        k = k+1;
        WEL{k} = [(1:NPER)' u*o(iw).LRC(i,:) o(iw).Q(:)*o(iw).fQ(i)];
    end
end

WEL = sortrows(cell2list(WEL));


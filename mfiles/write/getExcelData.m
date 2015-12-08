function [numHdr,numVals,txtHdr,txtVals,singleValue]=getExcelData(XLSname,SHTname,direction,varargin)
%GETEXCELDATA reads info from excel worksheet
%
% USAGE:
%    [numHdr,numValues,txtHdr,txtValues,singleValue] = ...
%           getExcelData( XLSname, SHTname, direction [,varName] [, options] );
%
% What it does:
%   reads data from excel file sheets, splitting the the sheet data into
%   headers and values, separating numerical values (fist two output
%   arguments) from the textual values (output arguments 3 and 4). A single
%   column of values can also be requested if its header is given an an
%   input argument.
%   The data headers are considered to be the last line of text labels
%   above the column data or the last column before the row data.
%   Blank lines and columns are automatically removed.
%
%   XLSname  'name of excel workbook file
%   SHTname  'name of sheet within workbook file'
%   direction  either 'H[orizontal]' or 'V[ertical] indicates hdr labels
%      direction in the worksheet.
%   varName   optional column header to request the value for a single
%      parameter.
%   options:
%        getExcelData( ... 'fill',true|false ... )
%              to fill empty cells with last value in same column
%        getExcelData( ... 'skip', vectorOfLinesToBeSkipped ... )
%        getExcelData( ... 'nanHdr', cellArrayWithNanInHdrReplacements ...)
%              headers must be text not numerical or formulas as they
%              yield NaN's. However you can set the 'nan' option to replace
%              these NaN's with the strings in the
%              cellArrayWithNanInHdrReplacements. Replacements of NaNs in
%              the hdr will be done with strings picked in sequence from
%              this cell array using its last value if more NaNs need to be
%              replaced than strings provided.
%
%
% TO 081231 110804 120415 130204 130626 151202

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin<3,
    error('getExelData:nargin:insufficientInput',...
        ['getExcelData: need at least three input arguments as follows from its usage:\n',...
        '[numHdr,numVals,txtHdr,txtVals,singleValue]=getExcelData(XLSname,SHTname,direction,varName)\n']);
end

% Usefull if fill is on and NaN's are desired
NaNthreshold = 1e50;  % values for which abs(values)>NaNthreshold will be set to NaN

% option to fill or not fill open field in worksheet.
[skip   ,varargin] = getProp(varargin,'skip',[]);
[fill   ,varargin] = getProp(varargin,'fill',true); % default is fill=true
[nanHdr ,varargin] = getProp(varargin,'nan' ,[]);
[varName,varargin] = getNext(varargin,'char',[]);

singleValueRequest = nargout>4 && ~isempty(varName);

if ~isempty(varargin)
    msgId = 'getExcelData:unusedArguments';
    warning(msgId,'Unused arguments in %s',mfilename);
end

%% Input checking to prevent trouble further down

% Third argument either H or V?
if isempty(regexp(direction,'^[hHvV]','once'))
   error('%s: The label direction (3rd argument) must be of type char and start with ''H'' or ''V''',mfilename);
end

%% first and second inputs must be of class char
if ~ischar(XLSname)
    error('%s: XLSname <<%s>> must be valid name of a Microsost Excel Spreadsheet',...
        mfilename,XLSname);
end
if ~ischar(SHTname)
    error('%s: SHTname <<%s>> must be a valid worksheetNm in workbook <<%s>>',...
        mfilename,XLSname,SHTname);
end

%% Second arg must be a worksheet in workbook given in first argument
[STATUS,SHEETS] = xlsfinfo(XLSname);
if isempty(regexp(STATUS,'Microsoft','once'))
    error('%s: XLSname %s not a legal %s',mfilename,'Microsoft Excel Spreadsheet');
end
if ~ismember(SHTname,SHEETS)
    error(['%s: SHTname %s not in workbook %s\n'...
           'REMEDY: maybe the workbook was not saved as an\n',...
           'Excel 5.0/95 workbook\n',...
           'The latter is necessary on non-Windows PC''s because\n',...
           'the Matlabs xlsread can''t read more recent fileformats of\n',...
           'Excel unless of version 5.0/95.\n',...
           'This will change depending on the Mathworks.'],...
           mfilename,SHTname,XLSname);
end

%% Ok, continue
direction = upper(direction(1));

warning('off','all');

if ismac
    warning('off','MATLAB:xlsread:Mode')
    try
    [~,~,Raw]=xlsread(XLSname,SHTname,'','basic');
        
    catch ME
   %     [STATUS,SHEETS] = xlsfinfo(XLSname);
        error(['%s: %s\n\n',...
            'REMEDY: First make sure the spelling of workbook and sheetName is exact.\n',...
            '   Notice that this spelling is case sensitive.\n',...
            '   Next, notice that this may happen on the MAC or other non-PC operating systems',...
            '   lacking excel.com so that xlsread must be run in basic mode.\n',...
            '   However, the Matlab function xlsread when running in basic mode\n',...
            '   will only work properly if the xls worksheet was saved\n',...
            '   as a <<Microsoft Excel 5.0/95 workbook>>.\n',...
            '   (On MS-Windows this (Matlab)-limitation does not apply).',...
            '   Therefore, saving the your Excel file in this format on Mac may solve this problem.\n',...
            '   You may use the command\n',...
            '   [~,sheets] = xlsfinfo(''%s'')\n',...
            '   to get a list of the sheets in the workbook that xlsread sees.\n',...
            '   My experience is that xlsread is not bug free, sometimes it does not see\n',...
            '   sheets that are definetly there. I have solved this by moving sheets around\n',...
            '   in the excel workbook (moving the invisible one to the front for instance,\n',...
            '   or copying the contents to a new sheet, remvoving the old one\n',...
            '   and renaming the new one to the old one.'],...
            mfilename,ME.message,XLSname);
    end
else
    try
        [~,~,Raw]=xlsread(XLSname,SHTname);
    catch ME
        error('%s: %s\nREMEDY: Check the exact spelling of the sheetName, this is case sensitive\n',...
            mfilename,ME.message);
    end
end

% eliminate lines to be skipped
if ~isempty(skip)
    Raw(skip,:) = [];
end

% Remove all empty columns
nans = cellfun(@(a) isnan(a(1)),Raw);
Raw  = Raw(~all(nans,2),:);

% Remove all empty columns
nans = cellfun(@(a) isnan(a(1)),Raw);
Raw  = Raw(:,~all(nans,1));


if direction=='H'

    txtCol  = all(cellfun( @(a) isnan(a(1)), Raw) | cellfun(@ischar,Raw),1);
    numCol  = ~txtCol;

    DataRow = all(cellfun(@isnumeric,Raw(:,numCol)),2);

    if ~any(DataRow)
        error(['%s: There are no numerical data rows in worksheet <<%s>> of workbook <<%s>>.\n',...
               'REMEDY: Make sure there is at least one numerical row in your worksheet.'],...               
            mfilename,SHTname,XLSname);
    end

    lblRow  = find(~DataRow);
    lblRow  = lblRow(end);

    if any(cellfun(@(a) isnan(a(1)), Raw(lblRow,:)))
        I = find(cellfun(@(a) isnan(a(1)),Raw(lblRow,:)));
        if isempty(nanHdr)                
            error(['There is at least one label that is EMPTY, or NUMERICAL or a FORMULA in\n', ...
                   'column <<%d>> in your worksheet <<%s>> of workbook <<%s>>.\n',...
                   'REMEDY: Verify that you don''t use character formulas in the labels in Excel (Matlab can''t read their textual values)'],...
                    I(1),SHTname,XLSname);
        else
            Raw(lblRow,I) = [nanHdr repmat(nanHdr(end),[1,numel(I)-numel(nanHdr)])];
        end
    end

    numHdr    = Raw(lblRow(end),numCol);
    numVals    = cell2mat(Raw(DataRow,numCol));
    txtHdr  = Raw(lblRow(end),txtCol);
    txtVals     = Raw(DataRow    ,txtCol);
        
    %% Fill in all NaN cells as of row 2
    % TO 130612
    if fill
       if any(isnan(numVals(1,:)))
            error(['When fill is on (default), Empty or NaN values are not allowed\n',...
                   'in the first row under the labels,\n',...
                   'in workbook <<%s>> sheet <<%s>> labels <<%s>>.\n',...
                   'Alternatively use value >%g to set cells to NaN''s'],...
                    XLSname,SHTname,sprintfs(' %s',numHdr(isnan(numVals(1,:)))),NaNthreshold);
       end

        % fill empty numeric data lines if necessary
        for iL = 2:size(numVals,1)
            I = isnan(numVals(iL,:));
            if ~isempty(I)
                numVals(iL,I)=numVals(iL-1,I);
            end
        end

        % fill empty text lines if necessary
        if ~isempty(txtHdr)
            for iL = 2:size(txtVals,1)        
                J = cellfun(@(a) isnan(a(1)), txtVals(iL,:));
                if ~isempty(J)
                    txtVals(iL,J) = txtVals(iL-1,J);
                end        
            end
        end
    end

    % potentially insert NaN by using an extreme threshold
    numVals(abs(numVals)>NaNthreshold)=NaN;

    if singleValueRequest
        I = strmatchi(varName,numHdr);
        if I(1)
            singleValue = numVals(:,I(1));
        elseif strmatchi(varName,txtHdr)
            singleValue = txtVals(:,strmatchi(varName,txtHdr));
        else
            error('Can''t find label <<%s>> in table in workbook <<%s>> worksheet <<%s>>.',...
                   varName,XLSname,SHTname);
        end
    end
    
    % Cut off any trailing NaNs from txtVals
    lastNonNaN = find(~all(cellfun(@(a) isnan(a(1)) ,txtVals),2),1,'last');
    txtVals(lastNonNaN+1:end,:)=[];

else % direction == 'V'

    txtRow  = all(cellfun(@(a) isnan(a(1)),Raw) | cellfun(@ischar,Raw),2);
    numRow  = ~txtRow;
    if ~any(numRow)
            error('Make sure there is at least one numerical row in your data of workbook <<%s>> worksheet<<%s>>.',...
                   XLSname,SHTname);
    end

    % Determination of DataCol requires at least columns with at least one numeric value
    % one numeric row (after the label columns)
    DataCol = any(cellfun(@isnumeric,Raw(numRow,:)),1);
    lblCol  = ~DataCol; % columns with all labels (need the last)
    i=find(~lblCol,1,'first');
    if ~isempty(i) && i>1
        lblCol(i:end)=false;
    end

    if any(cellfun(@(a) isnan(a(1)),Raw(:,lblCol)))
        I = find(cellfun(@(a) isnan(a(1)) ,Raw(:,lblCol)));
        if isempty(nanHdr)
            error(['There is at least one numerical label in row <<%d>> in your workbook <<%s>>, worksheet <<%s>>.\n',...
                   'Verify that you don''t use character formulas in the labels, as Matlab can''t read their textual output'],...
                    I(1),XLSname,SHTname);
        else
            Raw(I,lblCol) = [nanHdr repmat(nanHdr(end),[numel(I)-numel(nanHdr),1])];
        end
    end

    txtHdr  = Raw(txtRow,find(lblCol,1,'last'));
    txtVals     = Raw(txtRow,DataCol    );
    numHdr    = Raw(numRow,find(lblCol,1,'last'));

    % remove any cell with text from the numeric range
    for i=1:numel(numRow)
        for j=find(DataCol)
            if ischar(Raw{i,j})
                Raw{i,j}=NaN;
            end
        end
    end
    % then assign this now truly numeric range to numVals
    numVals    = cell2mat(Raw(numRow,DataCol));
    
    
    % Cut off any trailing NaNs from txtVals
    lastNonNaN = find(~all(cellfun(@(a) isnan(a(1)) ,txtVals),1),1,'last');
    txtVals(:,lastNonNaN+1:end)=[];

    if singleValueRequest
        I =strmatchi(varName,numHdr);
        if I(1)
            singleValue = numVals(I(1),:);
            singleValue(isnan(singleValue))=[];
        elseif strmatchi(varName,txtHdr)
            singleValue = txtVals(strmatchi(varName,txtHdr),:);
            singleValue(isnan(singleValue))=[];
        else
            error('Can''t find label <<%s>> in table in workbook <<%s>> worksheet <<%s>>.',...
                   varName,XLSname,SHTname);
        end
    end
end


function [o,PNTSRC]= getCin(o,basename,varargin)
% well = well.getCin(basename [{sheetNm,prefix1[,prefix2[,prefix3,...]]}])
% well = well.getCin(basename sheetNm);
% well = well.getCin(basename,prefixColHdr);
%
% get injection concentrations from worksheet sheetNm in workbook basename.
%
% basename is the workbook name (no need for extension .xls or .xlsx)
% sheetname is the sheet in the workbook that has columns with
%    concentration values per stess period for all requested species.
%    prefix is the prefix of the column headers with the concentration (or temp)
% Notice:
%    There can only be one sheetNm, the default is 'PER' in the workbook basename
%    There can be a prefix for every species. The default is 'C'
%      if no prefix is given. Then the prefix for the species will be C1,
%      C2 etc. And the column for the different wells C1_1 C1_2 ... C2_1,
%      C2_2 ..., where the first number indicates the species, and the
%      second the well number (i.e. well(iw).nr not iw.
%      Superfluous species names and wells will be ignored.
%   Combine species and well header prefixes in a cell array like
%     {'PER' 'temp' 'chloride'}
%     The sheetNm will be the name of these that matches any of the sheet
%     names in the workbook and defaults to 'PER'. The non-matching workds
%     are thus the names of the species in the order that they are
%     specified in SEAWAT and MT3DMS (and, therefore, in PNTSRC in the SSM
%     file.
%  example calls
%     well = well.getCin(basename,'PER');
%     well = well.getCin(basename,'PER','C');
%     well = well.getCin(basename,'PER','C','Temp','Pressure')
%       The concentrations are expected to be in columns in worksheet 'PER'
%       under header 'C', 'Temp', 'Pressure' with the wellNr attached to it, like in
%         C1, C2, Temp1 Temp2, Pressure1, Pressure2.
%
%    if there is only one column header that matches the prefix then this
%    column is used to supply data to all wells.
%
% SEE ALSO: gridObj/setWell wellObj/WEl wellObj/PNTSRC mfSetWells
%
%   TO 110426 120103 120408 130307 130412
%
% Copyright 2009-2012 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%% Number of species in the model
if nargin<2, return; end

[fill,    varargin] = getProp(varargin,'fill',true);
[indexFld,varargin] = getProp(varargin,'index','');
index = ~isempty(indexFld);

[MT3nams,MT3vals] =getExcelData(basename,'MT3D','Vert');
NCOMP = MT3vals(strmatchi('NCOMP',MT3nams),1);

[o.NCOMP] = deal(NCOMP);

% get prefixes and turn into cell immediately if char
[prefixes,varargin] = getNext(varargin,{'cell','char'},[]);
[prefixes,varargin] = getNext(varargin,'cell',prefixes);
if ischar(prefixes), prefixes = {prefixes}; end

if isempty(prefixes)
        sheetNm  ='PER';
        prefixes = {'C'};        
else
    [STATUS,SHEETS] = xlsfinfo(basename);
    if isempty(STATUS)
        error('%s: <<%s>> Does not seem a valid Microsoft Excel Wokrkbook',mfilname,basename);
    else % we have SHEETS as a cell array of strings
        I = ismember(prefixes,SHEETS);
        if any(I)
            if numel(find(I))>1
                error('%s: prefixes must have only one sheet name, you have %d',mfilename,numel(find(I)));
            else
                sheetNm = prefixes{I};
                prefixes(I)=[];
            end
        else
            sheetNm = 'PER';
        end
        if isempty(prefixes)
            prefixes = {'C'};
        end
    end
end

if ~isempty(varargin)
    fprintf('%s: some input not used.\n',mfilename);
end

%% Process the prefixes and get the data
[pernams,pervals,NPER]=getPeriods(basename,sheetNm,'fill',fill);

if NCOMP==numel(prefixes)
    % ok
elseif numel(prefixes)==1
    % ok, we use one prefix ofr all species   % strcmpi('C',prefixes)
else
    error(['%s: NCOMP=%d (see worksheet MT3D) does not match the number of\n',...
        'prefixes in the call <<%d>>.',...
        'If the number of species is greater than 1, specify prefixes of species explicitly in the call of wellObj.\n',...
        '(See help  wellObj)',...
        'Prefixes are strings naming the species to look for among the column headers of\n',...
        'your worksheet <<%s>> in workbook <<%s>>. The format of the header in must be\n',...
        'prefixN where prefix is the name of the species and N the well number.\n',...
        'For example: if you have two wells, nr 5 and 8 and the species salinity and temp,',...
        'you must have column headers in the PER worksheet:\',...
        '''salinity5'',''temp5'',''salinity8'',''temp8'' to specify salinity and temp values\n',...
        'for your wells 5 and 8 for all stress periods'],...
        mfilename,NCOMP,numel(prefixes),sheetNm,basename);
end

%% Check if prefix columns exist
% assert existance of columns with prefix of this component/species
for iComp = 1:NCOMP
    
    iCompMin = min(numel(prefixes),iComp);
    
    CCOL = strmatchi(prefixes{iCompMin},pernams); % not exact!
    
    if ~CCOL(1)
        % given prefix does not match any of the column headers, this is an
        % error.
        error('%s: No columns in with header prefix <<%s>> in worksheet<<%d>> of workbook <<%s>>',...
            mfilename,prefixes{1},sheetNm,basename);
    
    else
        % If no prefix was set, 'C' is the default. Then check to see if
        % prefix == 'C' if so make sure to not include CRCH and CEVT
        if strcmp('C',prefixes{iCompMin})
            I1 = strmatchi('CRCH',pernams);
            I2 = strmatchi('CEVT',pernams);
            if I1(1), CCOL=CCOL(~ismember(CCOL,I1)); end
            if I2(1), CCOL=CCOL(~ismember(CCOL,I2)); end
            if isempty(CCOL)
                error(['%s: There are no columns in with headers starting with prefix ''%s''\n',...
                    'for component nr <<%d>> in worksheet <<%s>> of workbook <<%s>>'],...
                    mfilename,prefixes{1},iComp,sheetNm,basename);
            end
        end
    end
    
    % At this point we have prefix(es).
    % Contiue processing each of the wells in turn.
    for iw=length(o):-1:1
        o(iw).NCOMP = NCOMP;
        if numel(CCOL)==1
            % all wells (1 or more) get their data from the same column
            cCol = CCOL;
        else
            
            % if there are more matching columns, look for columns for each
            % well in turn.
            
            try
                % see if there is a column matching the prefix together
                % with the well number, that is, with well(iw).nr.
                if ~index
                    cstr = sprintf('%s%d',prefixes{iComp},o(iw).nr);
                else
                    cstr = sprintf('%s%d',prefixes{iComp}.o(iw).UserData.(indexFld));
                end
                cCol = strmatchi(cstr,pernams,'exact');
                
                % if not, try same with underscore in name
                if ~cCol
                    if ~index
                        cstr2 = sprintf('%s_%d',prefixes{iComp},o(iw).nr);
                    else
                        cstr2 = sprintf('%s_%d',prefixes{iComp},o(iw).UserData.(indexFld));
                    end
                    cCol = strmatchi(cstr2,pernams,'exact');
                    if ~cCol
                        error('%s: Can''t find columns with prefix <<%s>> or <<%s>> in worksheet<<%s>> of workbook<<%s>>',...
                            mfilename,cstr,cstr2,basename,sheetNm);
                    end
                end
                
            catch ME
                % if not found, error
                error(['%s\n',...
                    'Can''t find column <<%s>> in worksheet <<%s>> of workbook <<%s>>\n',...
                    'corresponding to species <<%s>> for well(%d)'],...
                    ME.message,cstr,sheetNm,basename,prefixes{iComp},o(iw).nr);
            end
        end
        
        % found column for this well, so add concentration of iComp to this well
        o(iw).C(iComp,1:NPER)  = pervals(:,cCol)'; % as a horizontal vector for easy inspection
        o(iw).UserData.cCol(iComp) = cCol;
        o(iw).species{iComp} = prefixes{iCompMin};
    end
end
            

if nargout<2, return; end

%% Nargout>1, so that the PNTSRC array is requested
%  this array has [SP L R C Q] fields

% Count how many cells we have in total over all wells
% This requires that the wells have been set by o.toGrid(gr,HK)
% i.e. every o must have its grid information on board.

if isempty(o(iw).LRC)
    error('mfLab:wellObj_getQ:noLRCdata',...
    '%s: the wells have no grid info on board, run o=o.toGrid(gr,HK) first.',mfile);
end
    
LRC = vertcat([o.LRC]);
NCell = size(LRC,1);
n     = Inf;  % max number of cells per well to include in PNTSRC

%% All wells carry their grid  info and fow info. So we can proceed with PNTSRC

PNTSRC = cell(NCell,1);
u   = ones(NPER,1);

k=0;
for iw=1:length(o)
    Ncell = min(n,length(o(iw).idx));
    for icell = 1:Ncell
        k = k+1;
        if NCOMP==1
            PNTSRC{k} = [(1:NPER)' u*o(iw).LRC(icell,:) pervals(:,o(iw).UserData.cCol(1)) u*o(iw).ITYPE];
        else
            PNTSRC{k} = [(1:NPER)' u*o(iw).LRC(icell,:) pervals(:,o(iw).UserData.cCol(1)) u*o(iw).ITYPE pervals(:,o(iw).UserData.cCol)];
        end
    end
end

PNTSRC = sortrows(cell2list(PNTSRC));


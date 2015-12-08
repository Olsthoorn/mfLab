function [BCN,PNTSRC]=bcnZone(varargin)
%BCNZONE yields list data for stresses in MODFLOW and MT3DMS/SEAWAT
% bcnZone is a powerfull function to generate list input for MODFLOW
% stresses WEL, DRN, GHB, RIV, CHD and MT3DMS point sources PNTSRC as
% required by MT3DMS's SSM (source-sink mixing) package.
% bcnZone packs the lists in a cell array with one list per stress period.
% Each list containes line with
%  [SPnr L R C values]
%
% Notice: wellObj, MNWObj, MNW1Obj, MNW2Obj do not require such input
% lists, because these are generated automatically in the writeBCN and writeMNW
% functions directly from these objects.
%
% Example:
%      BCN         = bcnZone(basename,type,zoneArray,zoneVals,'fill',true|false)
%     [BCN,PNTSRC] = bcnZone(basename,type,zoneArray,zoneVals,concVals)
%     [BCN,PNTSRC] = bcnZone(basename,type,zoneArray,zoneVals,concVals,SP)
%     [BCN,PNTSRC] = bcnZone(basename,type,zoneArray,'SP',vectorOfActiveStressPeriods);
%
% BCN means "Boundary Conditions" (equivalent to stresses i.e.: WEL, DRN, GHB, RIV, CHD)
%     these are required by MODFLOW
% PNTSRC means "Point source". Required by the SSM package of MT3DMS/SEAWAT
%
% OUPUTS
%    BCN = array or array of cells containing the list input for a
%          particular boundary type (WEl, DRN etc).
%          The list has the following items on a line
%    [ stressPeriodNr  Layer Row Column values ]
%          Except for the stressPeriodNr, [Layer Row Column values] must be
%          exactly as required by the MODFLOW package of the type that is
%          specified.
%          The first argument, i.e. the stress period number, preceeds the
%          values required by the Modflow package in question. It makes
%          each input line/tuple unique, so that the order in which they are
%          presented to mfLab does not matter: mfLab will sort the input
%          accordingly and generate the input files for the different Modflow
%          packages in the right order.
%    PNTSRC = array similar to BCN, but required by MT3DMS/Seawat.
%          It contains the concentrations as follows
%    [ stressPeriodNr Layer Row Col Conc1 ITYPE Conc1 Conc2 Conc3 ...]
%          where ITYPE refers to a type of stress as defined in the MT3DMS
%          manual and Conc1 etc to the input conc of the corresponding
%          species for that cell and stress period. Note that input after
%          ITYPE is only required if the simulation is done with more than
%          one species.
%
% INPUTS
%    basename = the basename of the project used for (almost) all input and
%               output files
%    type     =is string indicating the stress type. It must be one of
%       'WEL' = well cells
%       'DRN' = drain cells
%       'DRT' = drain with return flow cells
%       'RIV' = river or stream cells
%       'GHB' = general-head boundary cells
%       'CHD' = (varying) constant heads
%       'MNW' = multi-node well cells
%       'MLS' = mass-loading source loading cells
%       'CCC' = constant-concentration cells
%       Notice: This list should be extended with each new stress type
%       being implemented.
%
%    zoneArray = an array of the size of the grid containing zone numbers.
%        zoneArray defines the spatial distribution of zones in the grid.
%
%    zoneVals = cell vector or cell array defining the input for the zones.
%        zoneVals has one row per zone as follows:
%        {zoneNr1 value value value ...;
%         zoneNr2 value value value ...
%         zoneNr3 value value value. }
%    zoneNri = the zoneNr partaining to the values that follow on the same line.
%        These values will be transferred to the input of the model and
%        must, therefore, match the inputs expected by the stress that is indicated
%        by the type argument.
%        For instance, WEL expects only the well flow. CHD expects two heads.
%        DRN requires an elevation and a conductance. See modflow the
%        manual for details.
%        regarding the way each value can be presented in the zoneVals line:
%           as a scalar: --> all cells with zoneNr get this value for all
%                   stress periods.
%           as a vector: with the same number of values as there are in the
%                   zone zoneNr in the model --> each cell gets the value
%                   pertaining to it, while the cells are the same for all
%                   stress periods.
%           as a string: the header in worksheet 'PER' in workbook [basename '.XLS]
%                   This yields one value for each stress period. This
%                   value will be the same for all cells in the zone.
%                   A different value for each zone and each
%                   stress period, requires a loop.
%           as a cell containing
%                   a scalar and a string, like {4 'stage'}
%                       In this case the different values in the column in
%                       worksheet'PER' under the header 'stage' will be
%                       multiplied by the scalar (here 4) in all cells
%                       pertaining to the zone in question.
%                   a vector and a string, where the length of the vector
%                     equals the number of cells in the zone pertaining to
%                     this value in zoneArray.
%                       In this case, the different values in the column in
%                       the worksheet 'PER' under the header given in the
%                       string will be multiplied by the values in the
%                       vector, which correspond exactly to the cells in
%                       the zonearray pertaining to this zone. The vector
%                       thus functions as a mutiplier array for the
%                       time-dependent values in the PER worksheet.
%
%     concVals: conc conc conc with one "value" per species.
%           conc can be a single value per species, or a vector with
%           the number of values equal to the number of cells in the
%           zone, or it is a string referring to the header of a column
%           in the PER sheet with one value per stress period.
%           In case we need a different value per cell and per stress
%           period, we have to use a loop.
%
%     SP or the string,value pair 'SP',vector can be used to specify stress
%     periods to be used instead of all stress periods. The number of
%     stress periods of the model is always used and obtained from the PER
%     worksheet.
%
%    'fill',true|false pair default true allows specifying whether
%     empty cells in the spreadsheet read will be filled from values above
%     or be set to NaN
%
% CHDDENSOPT = 2 % SEAWAT will use the environmental head (See seawat
% manual tm6A22.pdf from the SEAWAT download at http://water.usgs.gov/ogw/seawat/
%
% IMPORTANT: it is the user's responsibility to provide a
%    concentration for each component/species. mfLab cannot see of
%    you overrule the value NCOMP in sheet MT3D. It will,
%    therefore, just add any species concentrations to PNTSRC as
%    provided in concVals (=size(concVals,2)). This has no consequences
%    unless size(concVals,2) is less than the number of species
%    you are using, i.e. less than the actual NCOMP.
%
% See also: gridObj
%
% TO 120411 130814 (added spatially varying transient intput)

ITYPE = {'WEL','DRN','DRT','RIV','GHB','CHD','MNW','MLS','CCC'};


    if nargout>1 && nargin<5
        error('gridObj:bcn2:notEnoughInputArguments',...
            ['%s: to get PNTSRC you need 5 input arguments:\n',...
             'workbookName,bcnType,zoneArray,zoneVals,concVals\n',...
             'concVals missing\n'],mfilename);
    elseif nargin<4
        error('gridObj:bcn2:notEnoughInputArguments',...
       ['%s: insufficient input arguments concVals missing\n',...
        'Required are:\n',...
         ' 1:  workbookName with the PER worksheet to load data and number of stress periods.\n',...
         ' 2:  BCNtype, a string: ''WEL DRN DRT RIV GHB MNW MLR CCC''\n',...
         ' 3:  zoneArray like IBOUND with abs(zoneArray)==zone numbers\n',...
         ' 4:  zoneVals: one line per zone: [zoneNr zonvals]\n',...
         ' plus if nargout>1:\n',...
         ' 5:  concVals: one line per zone with one conc per compoment/species in MT3DMS or SEAWAT.\n',...
         ' zoneVals and or concVals may contain header names referring to corresponding columns in',...
         ' worksheet PER if individual data values for each stress period are required.\n',...
         ' in that case the repective array must be a cell array.'],...
         mfilename);
    end
    
%% assert input types

    [fill,varargin] = getProp(varargin,'fill',true);

    % get stress periods as an extra argument (undocumented)
    [SP,varargin] = getProp(varargin,'SP',[]);
    if isempty(SP)
        [SP, ~] = getNext(varargin,'double',[]);
    end
    
    [basename,varargin] = getNext(varargin,'char',[]);
    if ~ischar(basename)
        error('%s: First argument must be the basename of the workbook',mfilename);
    end
    
    [type,varargin] = getNext(varargin,'char',[]);
    if ~ischar(type)
        error('%s: Second argument must be string indicating the stress type like ''CHD'',''DRN'' etc.',mfilename);
    elseif ~ismember(type,ITYPE)
        error('%s: Second argument ''%s'' must be one of { %s }',mfilename,type,sprintfs(' ''%s''',ITYPE));
    end
    
    [zoneArray,varargin] = getNext(varargin,{'double','logical'},[]);
    if ~(isnumeric(zoneArray) || islogical(zoneArray))
        error('%s: Third argument must be a zoneArray like IBOUND',mfilename);
    end
    
    [zoneVals,varargin] = getNext(varargin,{'cell','double','logical','char'},[]);
    if isempty(zoneVals)
        error('Fourth argument must be zoneValues array of class cell or double (one line per zone)');
    end
    
    if nargout>1
        [concVals, varargin] = getNext(varargin,{'cell','double','char'},[]);
        if isempty(concVals)
            error(['%s: Fifth argument must be a conc array of class double or cell, one line per zone.\n',...
                'REMEDY: if arg is name of column in PER worksheet, make it a cell by putting it in { }'],mfilename);
        end
    end    

% stress or BCN types
    legalType={'CHD' 'WEL' 'DRN','RIV', 'GHB', 'MNW', 'MLS','CCC'};
    iSStype  =[  1     2     3     4      5      27    15    -1];
    if ~strmatchi(type,legalType)
        error('gridObj:bcn2:illegalITYPE',...
            ['%s: illegal type <<%s>>. Use one of:\n',...
            'CHD = constant-head cell',...
            'WEL = well\n',...
            'DRN = drain\n',...
            'RIV = river or stream\n',...
            'GHB = general-head boundary\n',...
            'MNW = multi-none well\n',...
            'MLS = mass-loading source\n',...
            'CCC = constant-concentration cell.'],mfilename,type);
    end
    
%% Algorithm
% We process one zone at a time, to which the corresponding line of the zoneValues cell array pertains
% nZones is the number of lines in the cell array

%% make sure every item is a cell. This allows uniform processing of
% different data types. For instance zoneVals could look like
%   { 3 4 h 'stage' c1; ...
%     [1 4 5] 4  4   {c2 'resistance2'}; ...
%     8 3  'head'  'stage2', 'resistance' }
% This means that 3 'zones are defined, ...
%      the cells with value 3 in zoneArray
%      the cells with values 1, 4 or 5 in zoneArray
%      the cells with value 8 in zoneArray.
% In the first line the 4 is a scalar, the h can be a scalar, vector or
% string, the 'stage' is a string and the c1 can be a scalar, vector or
% string. The values {c2 'resistance2'} on line 2, means that the
% resistance2 values that are given for the stress periods in the worksheet
% PER will be multiplied by c2, which is a scalar or a vector equal to the
% number of cells pertaining to the cell numbers 1, 4 or 5 in the
% zoneArray. This allows using c2 in this case as a multiplier array for
% the stress period values in the column 'resistance2' in the PER
% worksheet.M
%
% Notice that the number of arguments is the same on each line of zoneVals and
% is equal to the number required by the package indicated by the argument "type".

% Make sure that zoneVals is a cell array
if ~iscell(zoneVals)
    zV = zoneVals;
    zoneVals =cell(size(zoneVals));
    for i=1:numel(zoneVals)
        zoneVals{i} = zV(i);
    end
    clear('zV');
end


%% split the zoneNrs (first column of zoneValues) from subsequent data columns
nZones   = size(zoneVals,1);
zoneNr   = zoneVals(:,1);       % cell vector with zone numbers
zoneVals = zoneVals(:,2:end);   % zone without zone numbers
nCol     = size(zoneVals,2);

%% Get combined number of cells per zone of each zoneVals input line    

for iz=nZones:-1:1
    Iz{iz,1}      = find(ismember(zoneArray,zoneNr{iz}));
    Iz{iz}=Iz{iz}(:); % column vectors needed later
    if isempty(Iz{iz})
          error(['%s: zone(s) number <<%s >> was not found in the zoneArray!\n',...
              'Remedy: check your zone array.'],mfilename,sprintf(' %d',zoneNr{iz}));
    end
end
    
% The zoneVals can contain, scalars, a vector that is as
% long as there are cells with the zone in the zoneArray, or a string,
% which is the header in the PER worksheet, to select transient data, or,
% finally, zoneArray can be a cell containing a scalar or a vector combined
% with a string such as {4 'stage'}. In the latter case, transient values will be obtained from
% the corresponding column in the PER sheet under the header given in the
% string. And these values will be multiplied by the scalar or the vector
% given. This allows transient values to be given to the cells of the zone
% as a multiplier array times the transient value pertaining to a stress
% period.
    
%% check that all zoneVals only has scalars or strings or vectors of the length
% equal to the number of cells in each zone.
    for iz=1:nZones
        for i=1:nCol
            if ischar(zoneVals{iz,i}),                 continue; end
            if isscalar(zoneVals{iz,i}),               continue; end
            if isnumeric(zoneVals{iz,i})
                if sameSize(zoneVals{iz,i},zoneArray)
                    zoneVals{iz,i} = zoneVals{iz,i}(Iz{iz});
                end
            end
            if isvector(zoneVals{iz,i}) && isnumeric(zoneVals{iz,i}) && ...
                  ~(numel(zoneVals{iz,i}) == numel(Iz{iz}))
              error('vector zoneVals{%d,%d} must have length(Iz{%d}) = %d, not %d',...
                  iz,i,iz,numel(Iz{iz}),numel(zoneVals{iz,i}));
            end
            if iscell(zoneVals{iz,i})
                zoneVals{iz,i} = zoneVals{iz,i}(:)'; % make vector horizontal
                if numel(zoneVals{iz,i})~=2          % assert only two elements in vector
                    error('%s: cell zoneVals{%d,%d} must contain a string and a numeric value or vector',mfilename);
                end
                if ischar(zoneVals{iz,i}{1}), zoneVals{iz,i}=fliplr(zoneVals{iz,i}); end % put string last
                if ~ischar(zoneVals{iz,i}{end}), % assert string is present
                    error('%s: zoneVals{%d,%d} must contain a header of a column in the PER sheet',mfilename,iz,i);
                end
                % assert first element is scalar or vector of proper length
                if isnumeric(zoneVals{iz,i}{1})
                    % if entire 3D array is given select the cells
                    if sameSize(zoneVals{iz,i}{1},zoneArray)
                        zoneVals{iz,i}{1} = zoneVals{iz,i}{1}(Iz{iz});
                    end
                end
                if ~isscalar(zoneVals{iz,i}{1}) && ~(isvector(zoneVals{iz,i}{1}) && numel(zoneVals{iz,i}{1}) == numel(Iz{iz}))
                    error('%s: zoneVals{%d,%d} must contain a scalar or a vector of length Iz{%d} = %d, not %d',...
                        mfilename,iz,i,iz,numel(Iz{iz}),numel(zoneVals{iz,i}{1}));
                end
            end
        end
    end
   
%% Handle stress periods, SP, given as last (6th) argument of call, i.e.
%  when explicitly requested. Otherwise isempty(SP), in which case all
%  stress periods will be provided.

    [pernams,pervals,NT]=getPeriods(basename,'fill',fill);            

    if isempty(SP)
        SP = 1:NT;    % all stress periods
    else
        SP = SP(:)';  % explicityly specified stress periods.
        NT = numel(SP);   % number of explicitly given stress periods.
        
    % the number of stress periods cannot be zero.
        if isempty(SP) || any(diff(SP)<0)
            error('%s: Specified stress period numbers must increase.',mfilename);
        end
    end
    
%% allocate to store stress info
    % BCN will be a cell array of length NT. Each element contains a list
    % with length equal to the number of cells in the combined zones
    % each line of the list contains stress period number, L R C and the
    % values specified in zoneVals,
    
    % pointer index to before start and on end of the rows for each of the
    % zones in BCN and PNTSRC
    zoneEnd   = cumsum(cellfun(@numel,Iz));
    zoneStart = 1 + [0; zoneEnd(1:end-1)];

    % nCol is the number of specified values in zoneVals
    % zoneEnd(end) is total number of cells in all zones combined
    
    % Fill the last BCN cell as a template
    BCN{NT,1} = NaN(zoneEnd(end),4+nCol);
    for iz=nZones:-1:1
        BCN{NT}(zoneStart(iz):zoneEnd(iz),2:4) = cellIndices(Iz{iz},size(zoneArray),'LRC');
    end
    
    % Fill the last BCN cell first column with its stress period number    
    BCN{NT}(:,1)=SP(NT);  % add SP number in first column
    
    % Copy the last cell to all previous cells and replace the SP number
    % with the one pertaining to the cell to which the copying takes place
    for it=NT-1:-1:1
        BCN{it}     =BCN{NT};
        BCN{it}(:,1)=SP(it);  % add SP number in first column
    end
    
    
% At this point we have BCN withoutout the actual data, just the stress
% period number and the LRC of the cells have been filled out. The next
% step is to fill out the actual data, i.e. the values in zonVals and if
% applicable the concentrations.

%% Processing BCN first, process PNTSRC thereafter, if nargout>1        
    
%% Get the data for each zone and each column in zoneValues.
% This is transient data in principle.
    
    for iz=nZones:-1:1
        for j=nCol:-1:1 % for each of the values given for the zone

        % the zone values
        % v can be a scalar a vector or a string
        % this was asserted above.
            v=zoneVals{iz,j}(:);

        % expand 
        % if numeric, then transient values are all the same
            if isnumeric(v)
                
            % we may have a single value or a vector with one value for
            % each cell. In any case the data are then steady state.
                Values= ( v(:).*ones(size(Iz{iz})) ) * ones(1,NT);  % N*Nt should also work if v has a value for each zone
            elseif ischar(v)
                v = v(:)';
                % if v is string, look up the  header
                iCol=strmatchi(v,pernams);
                if ~iCol
                    error('%s: Can''t find header <<%s>> in PER worksheet, check zoneVals{%d,%d}\n',...
                        mfilename,v,iz,j+1);
                end
                Values =  ones(size(Iz{iz})) * (pervals(SP,iCol)');
                
            elseif iscell(v)
                % if v is cell, then v{1} is scalar or vector and v{2] a string/header
                iCol=strmatchi(v{2},pernams); % lookup header in PER sheet
                if ~iCol
                    error('%s: Can''t find header <<%s>> in PER worksheet, check zoneVals{%d,%d}\n',...
                        mfilename,v{2},iz,j+1);
                end
                % multiply numeric cell values by transient SP values
                Values= ( v{1}(:).*ones(size(Iz{iz})) ) * (pervals(SP,iCol)');
            else
                error('%s',mfilname);
            end

            % loop over the SP to fill the data into BCN
            % Here we have one corresponding value for each stress period for this zone
            for it=NT:-1:1
                BCN{it}(zoneStart(iz):zoneEnd(iz),j+4)=Values(:,it);
            end
        end
    end
    
    if nargout<2, return; end
    
%% Process all PNTSRC if narout>1
    
% convert concVals to cell array if necessary
    if ~iscell(concVals)  % could be all numbers, not needing PER sheet
        concVals = {concVals};
    end
    
% if necessary expand the size of the array to that of the zoneVals
    NL = size(concVals,1);
    for iz = nZones:-1:NL+1
        for j=1:size(concVals,2)
            concVals{iz,j} = concVals{NL,j};
        end
    end
    
%% Template P
    
% required size of P, the length of B is already ok
    nConc = size(concVals,2);
    if nConc==1
        nCol = 6;
    else
        nCol = 6+nConc;
    end

%% Set up the PNTSRC cell arrays
    PNTSRC{NT,1}         = NaN(zoneEnd(end),nCol);
    for iz=nZones:-1:1        
        PNTSRC{NT,1}(:,1:4)  = BCN{NT}(:,1:4);
    end
    PNTSRC{NT,1}(:,6)    = iSStype(strmatchi(type,legalType));  % column 6 == ITYPE
    PNTSRC{NT,1}(:,1)    = SP(NT);  % add stress period number in col 1
    for it=NT-1:-1:1
        PNTSRC{it,1}     = PNTSRC{NT,1};            
        PNTSRC{it,1}(:,1)= SP(it);  % add stress period number in col 1
    end
    
%% fill the PNTSRC, also by looping backward
    
    for iz=1:nZones
    % Process concentrations                  
        for j=nConc:-1:1 % {NrOfZones NrOfSpecies}

        % v may be string (PER column header)
        % so we must process each one in turn
        
        v=concVals{iz,j};
           
        % if numeric the data are steady state
            if isnumeric(v)
                Values= (v(:).*ones(size(Iz{iz}))) * ones(1,NT); % also works if v is vector
            
        % if not, a string, it refers to a header in the PER worksheet.
            else
                iCol=strmatchi(v,pernams,'exact');
                
            % no header is an error
                if ~iCol
                    error('gridObj:bcn2:wrongConcColHeader',...
                        ['gridObj/bcn2: Can''t find conc column header <<%s>> in PER worksheet\n',...
                         'check concVals{%d,%d array.'],v,iz,j);
                end
                
                Values= ones(numel(Iz{iz}),1) * pervals(SP,iCol)';
                if any(isnan(Values(:)))
                    warning('gridObj:bcn2:NaNinConc',...
                        ['gridObj/bcn2: Conc values in PER column <<%s>> has NaNs\n',...
                         'check if any cells in this column has no values (are empty).'],v);
                end
            end
            
            for it=NT:-1:1
                if j==1
                    PNTSRC{it,1}(zoneStart(iz):zoneEnd(iz),  5)=Values(:,it);  % c1 also goes in column 5 MT3D compatibility
                end
                
                if nConc>1
                    PNTSRC{it,1}(zoneStart(iz):zoneEnd(iz),j+6)=Values(:,it);
                end
            end
        end
    end
end
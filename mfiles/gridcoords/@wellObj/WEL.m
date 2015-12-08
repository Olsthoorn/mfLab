function WEL=WEL(o)
    % WEL = well.WEL() --- generate a WEL array for MODFLOW
    % where WEL is a cell array with one cell per stress period.
    % WEL{iPer}  = [iPer L R C Q]   -- if n is omitted
    % WEL{iPer}  = [iPer L R C Q wellNr] if nargin>1
    %
    % well = wellObj or an array of well objects.
    % The wells are assumed to carry the grid information (see wellObj.well)
    % and to have a vector of Q(1,1:NPER) on board.
    % NPER is obtained from the first well.
    %
    % TO 120512

    if isempty(o)
        warning('wellObj:WEL:empty',...
            'wellObj/WEL: input empty array of wellObj');
        WEL={};
        return;
    end
    
    %% Step 1, count how many well cells we have.
    % We generally have more than one well cell per well screen. Therefore,
    % we must count the total number of well cells before allocating
    % memory.
    
    N=0; % Number of well cells
    for iw=1:numel(o)
        N=N+size(o(iw).LRC,1);
        if isempty(o(iw).Q)
            error(['%s: Can''t set flow Q, because well(%d).Q==[].\n',...
                'REMEDY1: call wellObj with name of sheet with Q as 5th argument.\n',...
                'REMEDY2: call wellObj with 5th argument as cell array to specifying\n',...
                '   a) the sheetNm holding the stress period data with  Q, and\n',...
                '   b) the header of the column holding the Q values per stress period.\n',...
                'E.g. like this\n',...
                '   well = wellObj(basename,sheetNmWells,grid,HK,...\n',...
                '       {sheetNmHoldingStressPeriod ColNameOfQ[withOrWithoutWellNr]});\n',...
                '   well = wellObj(basename,''wells'',gr,HK,{''PER'',''Q_''});\n',...
                '   Notice that the header column name must match fprefix+wellNr\n',...
                '   or there is exactly one header that matches all the wells such as\n',...
                '   ''Q'' for Q1, Q2 Q3 ... or Q_ for Q_1 Q_2 ... or MyWell for Myell1 MyWell2, etc.'],...
                mfilename,iw);
        end
    end
    
    %% All wells must have a Q-vector of length NPER
    % Therefore, we can retrieve NPER from the first well
    NPER = numel(o(1).Q);
    
    %% Step 2: Allocate memory to hold LRC and Q for all these cells
    LRC         = NaN(N,3);    % Store LRC
    Q           = NaN(N,NPER); % Store flow for cells and stress period
    wellNr      = NaN(N,1);    % Well Nr

    %% Step 3: Populate these arrays
    k=0;
    for iw=1:numel(o)
       if size(o(iw).Q,2) ~= NPER
           error('wellObj:WEL:NoNPERQvalues',...
               'well %d (well Id=%d) does not have NPER Q values.',iw,o(iw).id);
       end
       
       m = size(o(iw).LRC,1);
       LRC(   k+(1:m),:)=o(iw).LRC;
       Q(     k+(1:m),:)=o(iw).fQ(:) * o(iw).Q;
       wellNr(k+(1:m),1)=o(iw).nr;
       k=k+m;
    end
    
    %% Step 4: Generate WEL cells
    
    for iPer=NPER:-1:1
        %% If nargin>1 used to signal that well Nr will be added to each line of WEL
        if nargin>1,
            WEL{iPer} = [ ones(N,1)*iPer LRC Q(:,iPer) wellNr];
        else
            WEL{iPer} = [ ones(N,1)*iPer LRC Q(:,iPer)];
        end
        WEL{iPer} = WEL{iPer}(~isnan(Q(:,iPer)) & Q(:,iPer)~=0, :);
    end
    
    WEL=WEL(~cellfun('isempty',WEL));
    

function PNTSRC=PNTSRC(o)
    % [WEL, PNTSRC] = well.PNTSRC() --- generate a PNTSRC array for MODFLOW
    % [ ~ , PNTSRC] = well.PNTSRC()
    %    PNTSRC     = well.PNTSRC()
    %
    % where PNTSRC is a cell array with one cell per stress period.
    % PNTSRC{iPer}  = [iPer L R C CSS ITYPE ISCC(1..NCOMP)]
    % n is optional number of columns if  n>6 (default)
    %
    % well = wellObj or an array of well objects.
    % The wells are assumed to carry the grid information (see wellObj.well)
    % and to have a vector of Q(1,1:NPER) on board.
    % NPER is obtained from the first well.
    %
    % recirculation wells are already in well.C and MNW1.C and MNW2.C
    % this functions just spits out the PNTSRC given the wells and
    % multi-node wells.
    %
    % TO 120512 121126
    
    % Step 1, count how many well cells to include.
    % That is numel(o(iw).idx) if we have wellObj and 1 if we have MNW1Obj or MNW2Obj
    
    N=0;  % Number of well cells
    for iw=1:numel(o)
        if o(iw).ITYPE==2  % i.e. class(o) == 'wellObj'
            N=N+numel(o(iw).idx);
        else % in the case of MNW only one record per well and SP
            N=N+1;
        end
    end
    
    %% All wells should have a Q vector of length NPER
    NPER = numel(o(1).Q);
    NCOMP= size(  o(1).C,1);
    
    if ~exist('n','var') || isempty(n) || n<6, n=6; end %#ok<NODEF>
    if NCOMP>1, n=6+NCOMP; end
    
    %% Step 2: Allocate memory to hold LRC and Q for all these cells
    PNTSRC{NPER,1}=NaN(N,n);  % One cell array per stress period
        
    LRC = NaN(N,3);      % Store LRC
    Q   = NaN(N,NPER);
    C   = NaN(N,NPER,NCOMP);   % Store flow for cells and stress period
    
    %% Step 3: Populate these arrays
    k=0;
    for iw=1:numel(o)
       if size(o(iw).C,2) ~= NPER
           error('wellObj:PNTSRC:NoNPERQvalues',...
               'well %d (well Id=%d) does not have NPER C values.',iw,o(iw).id);
       end
       
       if o(iw).ITYPE==2
           m = size(o(iw).LRC,1);
           LRC(k+(1:m),:) = o(iw).LRC;
           Q(  k+(1:m),:) = o(iw).fQ(:) * o(iw).Q;
           for icomp=1:NCOMP
                C(  k+(1:m),:,icomp)=ones(size(o(iw).fQ(:))) * o(iw).C(icomp,:);
           end
      else
           m = 1;
           LRC(k+(1:m),:) = o(iw).LRC(1,:);
           Q(  k+(1:m),:) = o(iw).Q;
           for icomp = 1:NCOMP
               C( k+(1:m),:,icomp) = o(iw).C(icomp,:);
           end           
       end
       k=k+m;
    end
    
    %% Step 4: Generate PNTSRC cells
    for iPer=1:NPER
        if NCOMP==1
            PNTSRC{iPer} = [ ones(N,1)*iPer LRC C(:,iPer,1) ones(N,1)*o(iw).ITYPE ];
        else
            PNTSRC{iPer} = [ ones(N,1)*iPer LRC C(:,iPer,1) ones(N,1)*o(iw).ITYPE YS(C(:,iPer,:))' ];
        end
        PNTSRC{iPer} = PNTSRC{iPer}(~isnan(Q(:,iPer)) & Q(:,iPer)~=0, :);
    end

    PNTSRC = PNTSRC(~cellfun('isempty',PNTSRC));
    
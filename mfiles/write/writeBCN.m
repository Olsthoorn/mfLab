function MXSS = writeBCN(basename,bcn,bcnp)
%WRITEBCN writes input file for MODFLOW's stress packages (WEL, DRN, RIV, GHB, CHD ...)
%
% Example:
%    MXSS  = writeBCN(basename,bcn,bcnp);
%
% Notice: bcn/BCN stands for "Boundary Condition" or "STRESS"
%
% The bcn.name must be one of 'WEL', 'MNW', 'GHB', 'DRN', 'RIV', 'CHD'
%
% MXSS is maximum number of sources and sink of the current boundary
%
% TO 070702 090807 100530 120410

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

%%

if ~isfield(bcn,'BCN')
    bcn.BCN = {};
elseif ~iscell(bcn.BCN)
    bcn.BCN = {bcn.BCN};
end

switch bcn.type
    case {'WEL','FLUX'}, minRecLen = 4; info={'Q'};
        %  if wellObj are present then first make a WEL list and add them to bcn struct
        %  this way wellObj are added to the list automatically and we don't have
        %  to worry about them
        if isfield(bcn,'well')
            if isempty(bcn.BCN)
                bcn.BCN = bcn.well.WEL';
            else
                bcn.BCN(end+(1:bcn.NPER)) = bcn.well.WEL';
            end
        end
    case 'DRN', minRecLen = 5; info={'Elevation', 'Conductance'};
    case 'RIV', minRecLen = 6; info={'Stage', 'Conductance', 'RivBotElev'};
    case 'GHB', minRecLen = 5; info={'Bhead', 'Conductance'};
    case 'CHD', minRecLen = 5; info={'Shead', 'Ehead'};
    otherwise,
        error('Illegal/unknown bcn.type=<<%s>>\n',bcn.type);
end

% line objects can be of any type
if isfield(bcn,'point')
    bcn.BCN = [bcn.BCN; bcn.point.BCN(bcn.gr,bcn.type,basename)];
end
if isfield(bcn,'line')
    bcn.BCN = [bcn.BCN ;bcn.line.BCN(bcn.gr,bcn.type,basename)];
end
% area objects can be of any type
if isfield(bcn,'area')
    bcn.BCN = [bcn.BCN; bcn.area.BCN(bcn.gr,bcn.type,basename)];
end


%% Are any data are provided for this package ?
if iscell(bcn.BCN)
    bcn.BCN=cell2list(bcn.BCN);
end

%% If river, then make sure stage is >= river bottom elevation
if strcmpi(bcn.type,'RIV')
    lgcl = bcn.BCN(:   ,5) < bcn.BCN(:   ,7);
           bcn.BCN(lgcl,5) = bcn.BCN(lgcl,7);
end

%% Guarantee that the head in RIV, GHB, DRN and CHD is above the bottom of the cell
%  in which the head is specified.
%
%  Cells with head below their bottom are moved to cell below.
%  This is not recursive, hence no perfec guarantee.
%  Alternative: eliminate these cells altogether or,
%  apply rule recursively and eliminate if below bottom of model.
%  TO 140311
%
if strmatchi(bcn.type,{'DRN','GHB','CHD','RIV'})
    
    I = bcn.BCN(:,1)>0; % skip asPrevious lines, indicacted by -iPer
    
    % Get elevation of bottom of cells in BCN list
    zB = bcn.gr.ZBlay(cellIndex(bcn.BCN(I,4),bcn.BCN(I,3),bcn.BCN(I,2),bcn.gr));

    % which lines have head or  stage < zB?
    K  = bcn.BCN(I,5)<zB;
    if sum(K)>0
        fprintf('Correcting %d %s cells with head below their bottom\n',sum(K),bcn.type);
    end
    
    % increase layer index to put the point in the next layer
    I=find(I);
    bcn.BCN(I(K),2) = bcn.BCN(I(K),2)+1;
end

%% Count the number of sources, which is needed when running MT3DMS or SEAWAT
MXSS = 0;

%% Determine the number of BCN/STRESS recores per stress period?

asPrevious=-1; % flag indicating to reuse BCN of previous stress period

% Number of input records per stress period
nDPer=zeros(bcn.NPER,1);
   
if ~isempty(bcn.BCN)
    for iper=1:bcn.NPER
        if any(bcn.BCN(:,1)==-iper) && ~any(bcn.BCN(:,1)==iper)
            % Instead of the number of records, record it as asPrevious
            nDPer(iper)=asPrevious;
            if iper==1,
                % asPrevious not allowed for first stress period.
                error('%s: asPrevious (negative stress period) cannot be used for first stress period !\n',mfilename);
            end
        else
            % Count the number of recoreds in this stress period
            nDPer(iper)=numel(find(bcn.BCN(:,1)==iper));
        end
    end
end

% Maximum number of active STRESS records at any one stress period
MXACTBCN=max(nDPer);
    
%% Write to MODFLOW input file for this stress
fid=fopen([basename, '.',bcn.ext],'wt');

%% 0
fprintf(fid,'# MATLAB writeBCN for %s,  %s\n',bcn.type,datestr(now));
fprintf(    '# MATLAB writeBCN for %s,  %s\n',bcn.type,datestr(now));

%% 1
if bcnp.NPar~=0
    fprintf(fid,'PARAMETER%10d%10d\n',[bcnp.NP,bcnp.MXL]);
end

%% 2
% printing the simulation header allowing for auxiliary variables
% while taking the exception for CHD into account

% ===== first print the variables =================================
fprintf(fid,'%10d',MXACTBCN);  % always

if ~strcmp(bcn.type,'CHD'), fprintf(fid,'%10d',bcn.ICB    ); end % all but CHD

%%
% Allow for AUXILIARY variables
if isfield(bcn,'AUX')
    Iaux = find(cellfun(@ischar,bcn.AUX));
    info = [info bcn.AUX(Iaux)];
    nAUX = numel(Iaux);
    for i=Iaux
        fprintf(fid,'  AUXILIARY %s',bcn.AUX{i});  % the variable names
    end
    
    %% if WEL, we use wellObj instead of a list. The list is generated from the
    %  wells using method wellObj/WEL above (line 31). Auxiliary variable
    %  values can, 'therefore', not be easily added to the end of the
    %  bcn.WEL records. Auxiiary variables for wells often comes down to
    %  IFACE, for which a good default is zero. This is assumed here. This
    %  is a limited approach, which should be improved in the future.
    %  We add nAUX columns of zeros to the bcn.BCN array:
    %  TO 130221
    if strcmpi(bcn.type,'WEL')
        bcn.BCN = [bcn.BCN zeros(size(bcn.BCN,1),nAUX)];
    end
       
else
    nAUX = 0;
end

% ===== next print information on these variables =================
fprintf(fid,'    <===== MXACTBCN');
if ~strcmp(bcn.type,'CHD'), fprintf(fid,'     I%sICB',bcn.type); end
if isfield(bcn,'AUX') && ~isempty(bcn.AUX),
    fprintf(fid,' +  Auxiliary variables');
end
fprintf(fid,'\n');
    
%% 3  always skipped because parameters are not (yet) implemented

if bcn.FREE
    fmt = ' %12.7g';
else
    fmt = ' %9.4g'; % works for numbers >-10e100, otherwise the size occupied will be 11 characters instead of 10 (=' %9g');
end

 for iPar=1:bcnp.NPar
    bcnp.NLST(iPar)=numel(find(bcnp.BCN(:,1)==iPar));
    fprintf(fid,'%11s%11s% 10d%10d\n',char(bcnp.PARNAM(iPar)),char(bcnp.PARTYP(iPar)),bcnp.Parval(iPar),bcnp.NLST(iPar));
    L=size(bcnp.parnams,2)-4;
    fprintf(fid,['%10d%10d%10d',repmat(fmt,[1,L]),'\n'],bcnp.BCN(bcnp.BCN(:,1)==iPar,2:end));
 end

%% The non parameter (ordinary) values for this boundary condition
if ~isempty(bcn.BCN)
    
    L=size(bcn.BCN,2)-4;

for iPer=1:bcn.NPER   %numel(bcn.ITMP)
    
   I2=find(bcn.BCN(:,1)==iPer);  % which BCN on in this period?

   if ~isempty(bcn.BCN) && iPer>1 &&  ...
           (nDPer(iPer) == asPrevious || ...
           (nDPer(iPer) == nDPer(iPer-1) && ...
                  all(all(bcn.BCN(I2,2:end)==bcn.BCN(I1,2:end)))))
                fprintf(fid,'%10d%10d     asPrevious {SP %d}\n',asPrevious,bcnp.NPar,iPer);
   else
       %5
       fprintf(fid,'%10d%10d     ITMP   NP  {SP %d}. Following lines: L R C ',...
           nDPer(iPer),bcnp.NPar,iPer,iPer); % ITMP NP  should be NP(iPer)
       for i = 1:numel(info), fprintf(fid',' %s',info{i}); end; fprintf(fid,'\n');
       
       MXSS = max(MXSS,nDPer(iPer)); % Max Nr of sources of this BCN
       %6
       %% Check record length
       if L+3<minRecLen + nAUX,
           error('%s: min rec length for %s is %d + NAUX = %d, total %d. You have %d',...
               mfilename,bcn.type,minRecLen,nAUX,minRecLen+nAUX,L+3);
       end
       
       if ~isempty(bcn.BCN) && nDPer(iPer)~=asPrevious,
            fprintf(fid,['%10d%10d%10d',...
               repmat(fmt,[1,L]),'\n'],bcn.BCN(I2,2:end)');
       else
           fprintf(fid,'%10d\n',asPrevious);
       end
   end
   I1=I2;
end

fclose(fid);
end
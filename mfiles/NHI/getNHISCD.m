function [well,surwell,kD,WEL]=getNHISCD(gr,fname,Ix,Iy)
%GETNHISCD retrieves the WEL stresses in NHI model on a per cel bases
%
% Example:
%    [well,surwell,kD,WEL] = getNHISLD -- NHI datafile for WEL
%
% Notice: the NHI model has no stess files for RIV, DRN, GHB. The data for
% these model stresses are given in areal ASCII files.
%
% INPUT:
%   gr     is the grid of the submodel, not the NHI model
%   fname  is the name of the SCD file of the NHI (I don't know what SCD means,
%          the file contains the wells of the NHI.
%   Ix, Iy are the NHI cell indices of the submodel we cut out of NHI
%
% The file has two types of wells: normal wells true coordinates and kD
% info specified on the same line, and sur wells without this information.
% I don't know wat sur means, but these wells seem te represent irrigation
% wells that exist areally in south east Brabant province. We output both
% well bundles separately as we do the kD struct and the subgrid.
%
% OUTPUT:
%   well     wellObj of normal wells
%   surwell  wellObj of sur or irrigation wells
%   kD       struct with kD at normal well sites
%   WEL      WEL{1} normal well list
%            WEL{2} sur well list
%
% TO 120429 120530

% This is the file holding all wells of the NHI
fprintf('Reading file    ''%s''\n',fname);
fprintf('The file holding all wells of the NHI.\n');

fid=fopen(fname,'r');  if fid<0, error('can''t open file <<%s>>',fname); end

% Skip first line
fgets(fid);

% Read number of lines (wells) in the file. Each line is a well, or rather
% a model cell representing an extraction by a well. Because the NHI uses
% one model layer for each aquifer, wells correspond practically always to
% a model layer, i.e. the screen penetrates a single model layer and,
% therefore, a single cell.
NwellsInNHI=fscanf(fid,'%d',1); fgetl(fid);

% This is where we are in the file (location of file pointer relative to
% begin of file.
p1=ftell(fid);

%% Set xGr and yGr of cell center
% The wells in the file are of two types.
%  --- wells with indiviual information connected to them
%  --- wells starting with "sur" and no further data. These wells are
%  probably irrigation wells and therefore merely produce an extra amount
%  of evapotranspiration locally. These sur wells are located mainly in
%  sandy areas in south-east Noord Brabant.
% Here we look up where these sur wells start. Then we can read in both
% series in the file separately so as to keep the detailed information of
% the first series, which the second series of wells lacks.

pCur=p1;  % file pointer of pointing at the current well

for iw=1:NwellsInNHI           % for all wells in the file
    s=fgets(fid);        % read line from file with well data
    pPrev=pCur;          % set pointer of previous well to current
    pCur =ftell(fid);    % update file pointer to current wel
    if strcmp(s(52:54),'sur');   % first irrigation well found ???
            % N is number of normal wells 
            N=iw-1;               % We need this to read them out hereafter)
            p2=pPrev;            % p2 is pointer to first sur well
        fseek(fid,0,'bof');      % set file pointer back to begin of file
        break;                   % jump out of for loop
    end
end

% we are now ready to start reading both series of well separately

%% Read first series of wells, those with extra info connected to them

%%
% move file pointer to just after the second line with the number of wells
fseek(fid,p1,'bof');

%%
% The following comment line shows the layout of the lines in the fill file
% up to the first sur well
%         5       319       671    -385.0        -1 ID =     10001; X =    167600; Y =    545450; kD =   1628.33;  code    5: bf and of in same wvp

%% Read N lines according to this format.
% The last expression is a regular expression telling to read whatever
% consiting of the characters in the list between te brackets (see regular
% expression in the doc of Matlab). This guarantees capturing the string up
% to the end of the input line.
W1 = textscan(fid,['%10d%10d%10d%10f%10d ID =%10d; X =%10f; Y =%10f; kD =%10f;  code%5d: ',...
    '%[abcdefghijklmnopqrstuvwxyz0123456789.,;!?:<=>" ]'],N);

%% Transfer the read-in data to actual wells

%% Get the wells of the first series (non sur wells)
% There are normal wells and sur wells in the file. I have no idea what sur
% stands for but these wells seem to be irrigation wells without a name,
% located area wide in south East Brabant.
% The first bundle of wells (non sur wells) have detailed information
% asscociated to them such as kD, x, y (see layout above)

%% Read out the first bundle of wells

% Layer Row Col indices of the first set of wells
LRC = [W1{1},W1{2},W1{3}];

LRCsub      = LRC; 
LRCsub(:,2) = LRCsub(:,2)-Iy(1)+1;
LRCsub(:,3) = LRCsub(:,3)-Ix(1)+1;

%%
% Get wells within the subgrid only
I = find(LRCsub(:,2)>0 & LRCsub(:,2)<=gr.Ny & LRCsub(:,3)>0 & LRCsub(:,3)<gr.Nx);

%%
% First allocate memory to store the wells
fprintf('Checking normal wells...\n');
if ~isempty(I)
    well(length(I),1) = wellObj;  % get NwellsInNHI well objects
    kD(  length(I),1) = kdObj;
end

%% Populate these wells and kD struct with the information from the file
%
% Allocate room for kD values in number equal to the number of wells (N) having
% the kD value associated to them
% Each record in the struct will hold the world x,y coordinates, the actual
% kD value and the Layer Row Col (LRC) of the cell, which is a well cell for
% this extra information is given in the file we are reading:
fprintf('Filling %d wells with data read from file:\n',NwellsInNHI);


%% Process for all wells of series 1 that are in the subgrid

[well.fQ]  = deal(1);           % fraction of total well Q from this cell
[well.Dt]  = deal(1);           % time step length
[well.t ]  = deal(2012);        % time of well (dummy)
[kD.remark]= deal('from NHI WEL file (april 2012)');

for i=1:length(I)    
    %%
    % Populate this well with the grid information
    % Outcommented itmes that cost too much time
    well(i).nr  = I(i);        % Use index in original well file
    kD(  i).nr  = I(i);
    
    well(i).x   = W1{7}(I(i));    % its world x-coordinate
    kD(  i).x   = well(i).x;
    
    well(i).y   = W1{8}(I(i));    % its world y-coordinate
    kD(  i).y   = well(i).y;
    
    well(i).ix  = LRCsub(I(i),3);      % its column index
    kD(  i).ix  = well(i).ix;
    
    well(i).iy  = LRCsub(I(i),2);      % its row    index
    kD(  i).iy  = well(i).iy;
    
    well(i).iLay= LRCsub(I(i),1);   % its layer index
    kD(  i).iLay= well(i).iLay;
    
    %well(i).idx = cellIndex( well(i).ix, well(i).iy, well(i).iLay, gr.size);
    %kD(  i).idx = well(i).idx;
    
    well(i).LRC = LRCsub(I(i),:);   % its Layer Row Col index vector
    kD(  i).LRC = well(i).LRC;
    
    well(i).Q   = W1{4}(I(i));    % well extraction (only steady state value is given in NHI data)
    
    well(i).remark = W1{end}{I(i)};  % Keep rest of input line as remark for this well

    kD(i).value = W1{ 9}(I(i));   % actual kD
    kD(i).code  = W1{10}(I(i));   % Don't know what this is
    
    %well(i).z   = [gr.ZTlay(well(i).idx), gr.ZBlay(well(i).idx)];
    %kD(  i).z   = well(i).z;
    
    %well(i).DZ  = abs(diff(well(i).z));
    %kD(  i).DZ  = well(i).DZ;
    
    %well(i).ztop= gr.ZTlay(well(i).iy,well(i).ix,1);
    
    %%
    % Show that wells are actually being processed
    if rem(i, 1000)==0, fprintf('.'); end
    if rem(i,25000)==0; fprintf('\n%d wells inspected, %d wells in wellbundle{1} in subgrid\n',length(W1{:}),length(I)); end
end
fprintf('\n%d wells inspected and %d wells in wellbundle{1} in subgrid\n',length(W1{1}),length(I));


%% Now get the sur wells and put them in wellbundle{2}
% sur wells are probably irrigation wells. The well file on the NHI.nu site
% do not specify what sur means.

% Set file pointer at first sur well
fseek(fid,p2,'bof');

%%
% Layout of the lines with sur wells
%         2       817      1064      6.25        -1 sur
% Reading till the end of the file with the layout of the sur wells
W2 = textscan(fid,'%10d%10d%10d%10f%10d%s','delimiter','\t');

%%
% Layer Row Col indices of the first set of wells
LRC = [W2{1},W2{2},W2{3}];

LRCsub      = LRC; 
LRCsub(:,2) = LRCsub(:,2)-Iy(1)+1;
LRCsub(:,3) = LRCsub(:,3)-Ix(1)+1;

%%
% Get wells within the subgrid only
Isur = find(LRCsub(:,2)>0 & LRCsub(:,2)<=gr.Ny & LRCsub(:,3)>0 & LRCsub(:,3)<gr.Nx);

%%
% First allocate memory to store the wells
fprintf('Checking sur wells ...\n');
if ~isempty(Isur)
    surwell(length(Isur),1) = wellObj;  % get NwellsInNHI well objects
end

for i=1:length(Isur)        
    %% Populate this well with grid data and its coordinates as before:
    surwell(i).nr=Isur(i);  % Use index in original well file
    surwell(i).LRC = LRCsub(Isur(i),:);
    surwell(i).x   = gr.xm(LRCsub(Isur(i),3)); 
    surwell(i).y   = gr.ym(LRCsub(Isur(i),2));
    surwell(i).ix  = LRCsub(Isur(i),3);
    surwell(i).iy  = LRCsub(Isur(i),2);
    surwell(i).iLay= LRCsub(Isur(i),1);
    surwell(i).Q   = W2{4}(Isur(i));
    
    % Global index for this well
    surwell(i).idx = cellIndex( surwell(i).ix, surwell(i).iy, surwell(i).iLay, gr.size);
    
    surwell(i).fQ  = 1;
    surwell(i).Dt  = 1;
    
    surwell(i).t   = 2012;
    
    surwell(i).created=now; % Stamp well with its creation time
    
    %%
    % Show that wells are being processed ...
    if rem(i, 1000)==0, fprintf('.'); end
    if rem(i,25000)==0; fprintf('\n%d wells inspected, %d wells found in wellbundle{2} in subgrid\n',length(Isur)); end
end    
fprintf('\n%d wells inspected and %d wells found in wellbundle{2} in subgrid\n',length(W2{1}),length(Isur));

fclose(fid);

%% Finalyze with summary
fprintf('\n%d wells inspected, %d wells found in subgrid\n\n',NwellsInNHI,...
    length(I)+length(Isur));

if ~exist(   'well','var') || isempty(I),       well=[]; end
if ~exist('surwell','var') || isempty(Isur), surwell=[]; end
if isempty(kD)   kD  =[]; end

WEL{2}=[];
if ~isempty(I),    WEL{1}=[ones(size(   well(:))) vertcat(   well.LRC) vertcat(   well.Q) ]; end
if ~isempty(Isur), WEL{2}=[ones(size(surwell(:))) vertcat(surwell.LRC) vertcat(surwell.Q) ]; end


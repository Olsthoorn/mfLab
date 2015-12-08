function [valuesTo,Ifr,Ito,dxTo] = gridsTransfer(xGrFr,valuesFr,xGrTo,code,dim)
%GRIDSTRANSFER % converts valuesFr corresponding to source grid, xGrFr, to valuesTo according
%   to toGrid, xGrTo.
%
%    works in any of the 3 spatial dimensions.
%
% Example:
%     valuesTo               = gridsTransfer(); % runs the self test
%    [valuesTo,Ifr,Ito,dxTo] = gridsTransfer(xGrFr,valuesFr,xGrTo,code)
%
%  valuesFr is te array in its original orientation
%  xGrFr and xGrTo are always in x-orientations, but their original orientation is
%  given in dim.
%
%  Hence:
%    numel(xGrFr) == size(valuesFr,dim)+1 and
%    numel(xGrTo) == size(valuesTo,dim)+1
%
%  code = 'geometric' implying geometric addition
%  or     'harmonic'  implying harmonic weighting
%  or     'abundant'  'median','max','min'
%  or     'divide',   'width' ,'thickness'  --- for cell widths (all mean the same)
%
% xGrFr, xGrTo are vectors that will be oriented along the x-axis (second dimension)
% no matther along which dimension they have been supplied to this function.
% Dim is the original orientaion of the values.
% The value array will be oriented along the x-axis before before being subjected
% to the algorithm.
% Before exiting the the function, the valuesTo will be oriented back along the original axis.
%
%  Algorithm:
%  The algorithm treats the data as being x-oriented; therefore, the
%  algorithm first permutes the xGrFr, xGrTo and valuesFr according to dim
%  along the x-axis.
%  The last step is back-permuation to the original dimension.
%
%  A intermittent grid, I, with values Vint is generated first.
%  Its cell boundaries, xGrInt, combines xGrFr and xGrTo so that
%  xGrFr and xGrTo correspond uniquelyinto xGrInt.
%  IFr and ITo hold indices from the intermediate grid into
%  the from and the to grid respectively. Their length equals the number of
%  cells in the intermediate grid and their values are the indices in
%  the respective grids.
%  The data from the from grid are then uniquely pulled into the
%  intermediate grid using
%     Vint(:) = VFr(IFr);  % no indexing required
%  Then the data from the intermittent grid must be transferred to the to grid.
%  Because each toGrid index may occur more than once in the ITo index vector, we
%  must do this tocell by tocell, with simultaneous weighting of the pulled values.
%  The weighting can be done geometrically, harmonically by median, max and
%  min, depending on the code passed to this function.
%
% Geometrically: ito is index over Vto
%     Vto(ito)= D(Ito==ito).*Vint(Ito==ito)/Dto(ito);  % note that DTo(ito)=sum(D(ITo==ito));
% Harmonically
%     Vto(ito)= 1./(D(Ito==ito)./Vint(Ito==ito)/Dto(ito))
% Abundant or median
%     Vto(ito)= median(Vint(Ito==ito),2);
% Minimum
%     Vto(ito)= min   (Vint(Ito==ito),[],2);
% Maximum
%     Vto(ito)= max   (Vint(Ito==ito),[],2);
%
%     target grid, thereby allowing either geometric or harmonic weighting.
%
% See also: modelObj
%
%  TO 120518
%  TO 130625 --> huge bug discovered and solved: DXint(J) instead
%                of DXint(Ifr(J)) in lines 207 & 209

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later


if nargin==0, selftest; return; end

permutes = {[2 1 3]; [1 2 3]; [1 3 2]};

%% Assert input

%% assert that dim is acceptable

if nargin<5 || ~ismember(dim,[1 2 3 'x' 'y' 'z'])
    error('You must specify dim and it must be one of[1 2 3 ''x'' ''y'' ''z'']');
end

% assert dim is one of [1 2 3]

if dim=='x', dim=2; elseif dim=='y', dim=1; elseif dim=='z', dim=3; end

%% assert xGrFr and xGrTo are vectors
if ~isvector(squeeze(xGrFr)), error('%s: xGrFr must be a vector'); end
if ~isvector(squeeze(xGrTo)), error('%s: xGrTo must be a vector'); end

% orient them along columns
xGrFr = xGrFr(:)';
xGrTo = xGrTo(:)';

% Prevent finding doubles by rounding to 4 decimals
xGrFr = roundn(xGrFr,4);
xGrTo = roundn(xGrTo,4);

% Remember if vectors are increasing or decreasing
sn = sign(diff(xGrFr([1 end])));

%% Assert xGrFr and xGrTo are either ascending or descending and unique
if ~all(sign(diff(xGrFr))==sign(diff(xGrFr([1 end]))))
    error('%s: xGrFr must be monotonically increasing or decreasing',mfilename);
end
if ~all(sign(diff(xGrTo))==sign(diff(xGrTo([1 end]))))
    error('%s: xGrTo must be monotonically increasing or decreasing',mfilename);
end
% Assert that both xGrFr and xGrTo are ascending or descending
if ~( sign(diff(xGrFr([1 end])))==sign(diff(xGrTo([1 end]))) )
    error('%s: xGrTo and xGrFr must be both increasing or decreasing',mfilename);
end

%% Assert that start and end of the inVecror and outVector match
if max(xGrTo)<=min(xGrFr)
    error(['%s: There is no overlap between the to and the from grids:\n',...
        'The max values of the toGrid [%g] must be > than the min value of the fromGrid [%g]\n',...
        'REMEDY: Check that the extends of your two grids match.'],...
        mfilename,max(xGrTo),min(xGrFr));
end
if min(xGrTo)>min(xGrFr)
    error(['%s: There is no overlap between the to and the from grids:\n',...
        'The min values of the toGrid [%g] must be greater or equal than the max value of the fromGrid [%g]\n',...
        'REMEDY: Check that the extends of your two grids match.'],...
        mfilename,min(xGrTo),max(xGrFr));
end

%% Assert that the dimension of valuesFr array in dim direction matches length of xGrFr
if size(valuesFr,dim) ~= numel(xGrFr)-1
    error('%s: size(valuesFr,dim) ~= numel(valuesFr)-1, i.e. %d~=%d-1',...
            mfilename,size(xGrFr,2),size(valuesFr,2));
end

%% Assert correct code for weighing the toValues
code =lower(code);

switch code
    case {'k','g'}, code='geometric';
    case {'c','h'}, code='harmonic';
end

validCodes = {'harmonic','geometric','abundant','median',...
              'mininum','maximum','divide','thickness','width'};
            
I = strmatchi(code,validCodes);
if ~I(1)
    error('%s: code must be one of %s\n',mfilename,sprintfs(validCodes,', %s'));
else
    code = validCodes{I(1)};
end

%% orientate valuesFr along the column axis, i.e. force dim to be true (dim==2)
valuesFr = permute(valuesFr, permutes{dim});

% use xdim below instead of dim and remember dim for back orientatio at the end
xdim=2;

%% Algorithm

%% Generate the intermediate grid xGrInt, which is a merge of xGrFr and xGrTo
%  to ensure all valuesFr will have their place within it.
xGrInt  = unique([xGrFr xGrTo]);

%% If the original grids were descending, then flip them to guarantee ascending vectors

if sn<0,
    xGrFr    = fliplr(xGrFr);
    xGrTo    = fliplr(xGrTo);
    valuesFr = valuesFr(:,end:-1:1,:);
end

dxTo  = abs(diff(xGrTo));
dxInt = abs(diff(xGrInt));
%% intermediate grid cell centers, will never coincide with cell boundaries
%  They are use to link the cells of the xGrTo and xGrFr

xmInt  = 0.5*(xGrInt(1:end-1)+xGrInt(2:end));  % intermediate cell centers

%% Pointers from intermediate grid into source grid Ifr and into target grid Ito
% Non-overlapping source-merged or merged-target grids get NaN as pointers
% Need to work in steps because min(NaN,x) --> x and not NaN !!
Ifr = floor(interp1(xGrFr,(1:length(xGrFr))',xmInt)); % index grid 1 in merged
Ifr(Ifr>length(xGrFr))=length(xGrFr);

Ito = floor(interp1(xGrTo,(1:length(xGrTo))',xmInt)); % index grid 2 in merged
Ito(Ito>length(xGrTo))=length(xGrTo);

% find start and end index of if to grid in intermediate grid
[~,First] = unique(Ito,'first');  % First where they start in vInt
[~,Last ] = unique(Ito,'last');   % and Last where they end

%% Transfer data from source grid to intermediate grid

DXint = repmat(dxInt,[size(valuesFr,1),1,size(valuesFr,3)]);

Vint  = NaN(size(DXint));     % Vint are the values in the intermedate grid
J     = find(~isnan(Ifr)); % J are valid indices in Ifr in case of non overlapping grid portions

if ismember(code,{'divide','thickness','width'})
    % This is necessary for cell widths,
    DXfr = repmat(abs(diff(xGrFr)),[size(valuesFr,1),1,size(valuesFr,3)]);    
    Vint(:,J,:) = valuesFr(:,Ifr(J),:) .* DXint(:,J,:)./DXfr(:,Ifr(J),:);
    code = 'width';
    %% Check balance, total thickness must not change
   % [sum(Vint,2) sum(valuesFr,2)]
else
    DXto = repmat(abs(diff(xGrTo)),[size(valuesFr,1),1,size(valuesFr,3)]);    
    switch code
        case 'geometric'
                Vint(:,J,:) = bsxfun(@times,DXint(J),valuesFr(:,Ifr(J),:));  % include intermediate dx as weight
        case 'harmonic'
                Vint(:,J,:) = bsxfun(@rdivide,DXint(J),valuesFr(:,Ifr(J),:));  % include intermediate dx as weight
        otherwise
            %error('%s: unknown case %s',mfilename,code);
    end
end
%% Fill new grid
% Non-overlapping portion initially zero to allow accumulation
valuesTo=zeros(size(valuesFr,1),size(xGrTo,2)-1,size(valuesFr,3));

% We will pull the data from the intermedate grid into the target grid.
% This way we can deal with integer grid, by selecting the most abundant
% value of the integers in the intermediate grid where it overlaps cells of
% the target grid. Most abundant is median or round(median). This can also
% be applied to zone arrays. However some data will drop out like small
% rows with boundaries. Alternatives are choosing minium or maximum values
% in a range. This is up to the user.
% Non-overlapping is autmatically taken care of. It yieds NaN's in the
% resulting valuesTo except if min or max are used.
for j=1:size(valuesTo,xdim)  % J is overlapping part of source<->intermediate
    switch code
        case 'width'
            valuesTo(:,j,:) =  sum(Vint(:,First(j):Last(j),:),xdim);            
        case  'geometric',
            valuesTo(:,j,:) =  sum(Vint(:,First(j):Last(j),:),xdim)./DXto(j);
        case 'harmonic'
            valuesTo(:,j,:) =  DXto(j)./sum(Vint(:,First(j):Last(j),:),xdim); 
        case {'abundant','median'}
            valuesTo(:,j,:) =  round(median(median (Vint(:,First(j):Last(j),:),xdim),xdim));
        case 'maximum'
            valuesTo(:,j,:) =  max(Vint(:,First(j):Last(j),:),[],xdim);
        case 'minimum'
            valuesTo(:,j,:) =  min(Vint(:,First(j):Last(j),:),[],xdim);
    end
end

%% If the original array has been flipped make sure to reflip the new one

if sn<0,
    valuesTo = valuesTo(:,end:-1:1,:);
end

%% Finally valuesTo into the original orientation of the valuesFr
valuesTo = permute(valuesTo,permutes{dim});



function selftest

%% Check that total is preserved
xGrFr = unique(100*[0 rand(1,18) 1]);
xGrTo = unique(100*[0 rand(1, 7) 1]);
valuesFr = rand(10,size(xGrFr,2)-1,5);
DXfr = repmat(diff(xGrFr),[size(valuesFr,1),1,size(valuesFr,3)]);
DXto = repmat(diff(xGrTo),[size(valuesFr,1),1,size(valuesFr,3)]);

permutes={[2 1 3]; [1 2 3]; [1 3 2]};

%% To check that harmonic also works set all values of Fr equal to 1
valuesFr(:,:,:)=1;

%% To verify for other dimensions
dim = 3;

xGrFr    = permute(xGrFr   ,permutes{dim});
xGrTo    = permute(xGrTo   ,permutes{dim});
valuesFr = permute(valuesFr,permutes{dim});
DXfr     = permute(DXfr    ,permutes{dim});
DXto     = permute(DXto    ,permutes{dim});

%% Geometric
valuesToGeometric = gridsTransfer(xGrFr,valuesFr,xGrTo,'g',dim);

%% Harmonic
valuesToHarmonic  = gridsTransfer(xGrFr,valuesFr,xGrTo,'h',dim);

xGrFr    = permute(xGrFr   ,permutes{dim}); %#ok<NASGU>
xGrTo    = permute(xGrTo   ,permutes{dim}); %#ok<NASGU>
valuesFr = permute(valuesFr,permutes{dim});
DXfr     = permute(DXfr    ,permutes{dim});
DXto     = permute(DXto    ,permutes{dim});

%% Check or verify

display([sum(DXfr.*valuesFr,2) sum(DXto.*valuesToGeometric,2)]);
display([sum(DXfr.*valuesFr,2) sum(DXto.*valuesToHarmonic ,2)]);


%% Second self test

xGrFr = [1 2 3 4 5 6];
%xGrTo = [1 1.5 2 2.25 2.5 2.75 3 4  5.25 6];
xGrTo = 1 + unique([0 5*rand(1,20) 5]);
dx =diff(xGrFr);

dxTo = gridsTransfer(xGrFr,dx,xGrTo,'divide',2);
sum(dxTo)



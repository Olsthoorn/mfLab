function varTo = gridsTransfer(wfr,varFr,wto,code,dim)
%% varTo = gridsTransfer(wfr,varFr,wto,code,dim)
%  Allows to resample varFr defined on linear coordinates wfr onto a nuew
%  linear grid wto, yielding varTo in the same size and dimesion as varrFr
%  perpendicular to dim. Hence varFr may be a 3D block while wto and wfr
%  are vectors in dimension dim.
%
%  wto is subject to the following conditions:
%  *  the outer planes of both grids must coincide.
%  *  the two dimensions not tranferred must have the same number of elements
%  in the two grids.
%  * the intermediate planes of Znew have arbitrary elevations but must be in
%    sequence from high to low.
%
%  use Znew = wfr.newZ(layersToKeep,subdivisions) to generate Znew
%
%  code = 'geometric' implying geometric addition
%  or     'harmonic'  implying harmonic weighting
%  or     'abundant'  'median','max','min'
%  or     'divide',   'width' ,'thickness'  --- for cell widths (all mean the same)
%
%  Procedure:
%  all cells are put in a single vector in sequence along the dimension to
%  be tranferred.
%  A intermediate grid wint and vInt with values vInt are generated first.
%  Its cell boundaries, wInt, combine those of the two grids
%
%  The data from the grid are then pulled into the intermediate grid
%  using V(:) = VFr(IFr), with no indexing required
%
%  Then the data from the intermiate grid vInt must be transferred to the to output grid varTo.
%  Because each toGrid index may occur more than once in the output index vector, we
%  must do this cell by tocell with simultaneous weighting of the pulled values.
%  Weighting is specified by the 'code' argument and can be done geometrically,
%  harmonically by median, max and min.
%
% Geometrically: ito is index over Vto
%  Vto(ito)= D(Ito==ito).*V(Ito==ito)/Dto(ito);  % note that DTo(ito)=sum(D(ITo==ito));
% Harmonically
%  Vto(ito)= 1./(D(Ito==ito)./V(Ito==ito)/Dto(ito))
% Abundant or median
%  Vto(ito)= median(V(Ito==ito),2);
% Minimum
%  Vto(ito)= min   (V(Ito==ito),[],2);
% Maximum
%  Vto(ito)= max   (V(Ito==ito),[],2);
%
%  target grid, thereby allying either geometric or harmonic weighting.
%
% Copyright 2009-2012 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later
%
% See examples/tutorial/modelObj for verification and application
% See examples/swt_v04/SWIM22/ASTRseasonal for application
%
%  TO 120821

if any(wfr.LAYCBD>0), error('%s: LAYCBD must be all zero (no confining beds'); end

%% Algorithm
code =lower(code);

switch code
    case {'k','g'}, code='geometric';
    case {'c','h'}, code='harmonic';
end

%% Check size of varFr
 if ~strcmpi(code,'layer') && ~numel(varFr)==numel(wfr)-1
     error('%s: size of varFr [%d] does notmatch size of input vector [%d]',mfilename,numel(varFr),numel(wfr)-1);
 end

%%
% Combine all columns of the data and Z-arrays into column vectors and
% transfer it back after the data have been transferred to the new grid.
% This way we can handle the entire grid in a linear interpolation
% operation. Care must be taken that we remain consistent with our
% numbering. So it is not allowed to throw away output grid planes, but
% thowing away input grid duplicates is no problem but desired.
% We dont' care about the input numbering, only the output numbering
% matters.

%% To grid and to values

% elevation vector grid average
wfr = wfr(:);
wTo = wto(:);
wInt = unique([wfr; wTo]);
wmInt=0.5*(wInt(1:end-1)+wInt(2:end));
    
% indices from intermediate grid wInt into fr and to grid
Ifr = min(numel(wfr)-1,floor(interp1(wfr,1:length(wfr),wmInt))); % index grid 1 in merged, size(Nint,1)
Ito = min(numel(wTo)-1,floor(interp1(wTo,1:length(wTo),wmInt))); % index grid 2 in merged, size(Nint,1)

% find start and end index of if to grid in intermediate grid
[~,First] = unique(Ito,'first');  % First where they start in vInt
[~,Last ] = unique(Ito,'last');   % and Last where they end

if strcmpi(code,'zlist')
    % one value per layer, don't need to transfer entire grid, which may be large
    varTo = zeros(size(wto(1:end-1)));
    
    % pull values from from grid into intermediategrid
    vInt = varFr(Ifr);
    
    % pull values from intermediate grid into to-grid
    for j=1:length(varTo)
        varTo(j) = round(median(vInt(First(j):Last(j))));
    end
    return;
end

%% working with entire grids
%% put the well size in a single vector after permuting the grids the grid
% to aline the transfer diretion with the primary dimension (dimension 1)
% along the first dimension.
dFr  = diff(wfr);
dTo  = diff(wto);
switch dim
    case {1,'y'}
        varFr=permute(varFr,[1,2,3]);
    case {2,'x'}
        varFr=permute(varFr,[2,1,3]);
    case {3,'z'}
        varFr=permute(varFr,[3,2,1]);
end
varFr=varFr(:);

varTo= zeros(size(dTo));

%% Intermediate grid and intermediate values
wInt  = unique([wfr;wTo]);                 % intermediate grid-line coords
wmInt = 0.5*(wInt(1:end-1)+wInt(2:end));   % intermediate center coords
dInt =  diff(wInt);                        % intermediate grid cell widths
vInt  = zeros(size(dInt));                 % intermediate grid values

%% Pointers from intermediate grid into from grid Ifr and into target grid Ito
Ifr = min(numel(wfr)-1,floor(interp1(wfr,1:length(wfr),wmInt))); % index grid 1 in merged, size(Nint,1)
Ito = min(numel(wTo)-1,floor(interp1(wTo,1:length(wTo),wmInt))); % index grid 2 in merged, size(Nint,1)

%% Get unique varTo indices
[I,First] = unique(Ito,'first');  % First where they start in vInt
[~,Last ] = unique(Ito,'last');   % and Last where they end

%% Use only the overlapping grids
J  = find(~isnan(Ifr)); % J are valid indices in Ifr in case of non overlapping grid portions

if strcmpi(code,'layer')
    %% In this case the we deal with a stress list with LRC indices of which
    % the index for the current dimension has to be adapted to the new grid.
    % We pass the linear layer index. We rcompute it for the permuted
    % grid. The index then correspoinds with the
    % linear cell in the wfr vector. We can thus immediately compute
    % the center of the cell of the stress in the wfr vector. This
    % w-value (grid coordinate) is the same in the zTo vector, which is in
    % sequence of the idx of the new grid. We kan thus immediately look up the
    % linear index of the zNew array, compute the new layer and
    % transfer that back to the stress list.
    % The condition for this to work is that the top and
    % bottom planes of the from and the to networks are the same.

    % varFr contains the linear index of the cells
    CRL = cellIndices(varFr,wfr.size,'CRL');
    switch dim
        case 1
            idx = cellIndex(CRL([1 2 3]),[wfr.Nx,wfr.Ny,wfr.Nz]); 
            ymIdx = 0.5*(wfr(idx)+wfr(idx+1));
            varTo = min(floor(interp1(wTo,1:length(wTo),ymIdx)),length(wTo)); % linear index grid 2, size(Nint,1)
        case 2
            idx = cellIndex(CRL([2 1 3]),[wfr.Ny,wfr.Nx,wfr.Nz]); 
            xmIdx = 0.5*(wfr(idx)+wfr(idx+1));
            varTo = min(floor(interp1(wTo,1:length(wTo),xmIdx)),length(wTo)); % linear index grid 2, size(Nint,1)
        case 3 % varFr = idxOld, varTo =idxNew
            idx = cellIndex(CRL([3 1 2]),[wfr.Nz,wfr.Nx,wfr.Ny]); 
            zmIdx = 0.5*(wfr(idx)+wfr(idx+1));
            varTo = min(floor(interp1(wTo,1:length(wTo),zmIdx)),length(wTo)); % linear index grid 2, size(Nint,1)
    end
    return;
end


%% Fill vInt (pull varFr values into vInt
if ismember(code,{'divide','thickness','width'})
    % This is necessary for cell widths,
    vInt(J) = varFr(Ifr(J)) .* dInt(J)./dFr(Ifr(J));
    code = 'width';    
else % just pull the value from varFr into vInt
    vInt(J) = varFr(Ifr(J));
end

%% Fill new to-grid
% We will pull the data from the intermedate grid into the target grid.
% This way we can deal with integer grid, by selecting the most abundant
% value of the integers in the intermediate grid where it overlaps cells of
% the target grid. Most abundant is median or round(median). This can also
% be applied to zone arrays. However some data will drop out like small
% rows with boundaries. Alternatives are choosing minium or maximum values
% in a range. This is up to the user.
% Non-overlapping is autmatically taken care of. It yieds NaN's in the
% resulting var2 except if min or max are used.

switch code
    case {'width','layer'}
        for j=I'
            varTo(j) = sum(vInt(First(j):Last(j)));
        end
    case  'geometric',
        vInt=dInt.*vInt;
         for j=I';
             varTo(j) = sum(vInt(First(j):Last(j)))./dTo(j);
         end
    case 'harmonic'
        vInt = dInt./vInt;
        for j=I'
            varTo(j) = dTo(j)./sum(vInt(First(j):Last(j)));
        end
    case {'abundant','median'}
        for j=I'
            varTo(j) = round(median(vInt(First(j):Last(j))));
        end
    case 'maximum'
        for j=I'
            varTo(j) = max(vInt(First(j):Last(j)));
        end
    case 'minimum'
        for j=I'
            varTo(j) = min(vInt(First(j):Last(j)));
        end
    case 'zlist'
        for j=I';
            varTo(j) = vInt(First(j));
        end
    otherwise
        legalCodes ={'harmonic','geometric','abundant','median',...
            'mininum','maximum','divide','thickness','width',...
            'layer','zlist'};
        if ~ismember(code,legalCodes)
            error(sprintf('%s: code must be one of\n%s\n',mfilename,repmat(' ''%s''',[1,length(legalCodes)])),legalCodes{:});
        end
end        

%% Round off by reshaping varTo to its  (Ny,Nx,Nz) form
switch dim
    case 1
        varTo = permute(reshape(varTo,[wto.Ny,wto.Nx,wto.Nz]),[1 2 3]);
    case 2
        varTo = permute(reshape(varTo,[wto.Nx,wto.Ny,wto.Nz]),[2 1 3]);
    case 3
        varTo = permute(reshape(varTo,[wto.Nz,wto.Nx,wto.Ny]),[3 2 1]);
end


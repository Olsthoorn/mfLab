function varTo = gridTransfer(grOld,varFr,grNew,code,dim)
%% var2 = gridObj.gridTransfer(var,grNew,code,dim)
%  Allows to resample a complete grid into a new grid implied by planes and grOld to grNew.
%  This is a powerfull, sophisticated and advanced procedure because it
%  works with arbitrary grids of the MODFLOW type as long as no confining beds are present.
%  Its use is to adapt grids so that they can be used with transport models
%  like mt3dms and seawat.
%
%  grNew is subject to the following conditions:
%  *  the outer planes of both grids must coincide.
%  *  the two dimensions not tranferred must have the same number of elements
%  in the two grids.
%  * the intermediate planes of Znew have arbitrary elevations but must be in
%    sequence from high to low.
%
%  use Znew = grOld.newZ(layersToKeep,subdivisions) to generate Znew
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

% factor to remove lines if less than MINDZ/10 apart to prevent ghost doubles when
% using interp1
ROUND = 10/grOld.MINDZ; 

if any(grOld.LAYCBD>0), error('%s: LAYCBD must be all zero (no confining beds'); end

%% Algorithm
code =lower(code);

switch code
    case {'k','g'}, code='geometric';
    case {'c','h'}, code='harmonic';
end

%% check dim
switch dim
    case {1,'y','Y'}, dim=1;
    case {2,'x','X'}, dim=2;
    case {3,'z','Z'}, dim=3;
    otherwise
        if isnumeric(dim), dim = sprintf('%d',dim); end
        error('%s: illegal dim value (%s), must be 1,2,3, ''x'',''y'',''z''',mfilename,dim);
end
%% Check size of varFr
 if strcmpi(code,'zlist') && ~(length(varFr)==grOld.Nlay)
     error('%s: length of varFr [%d] does not match grid.Nlay [%d]',mfilename,length(varFr),grOld.Nlay);
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

%% Check the vertical extent of the grid. Should perhaps be done over all columns
%  separately in case layers are not uniform

delta =1e-3;
switch dim
    case 1
        if abs(grOld.yGr(  1)-grNew.yGr(  1))>delta || ...
           abs(grOld.yGr(end)-grNew.yGr(end))>delta
            error('%s: yGr(1) and yGr(end) must be equal in both grids',...
                mfilename);
        end
    case 2
        if abs(grOld.xGr(  1)-grNew.xGr(  1))>delta || ...
           abs(grOld.xGr(end)-grNew.xGr(end)>delta)
            error('%s: xGr(1) and xGr(end) must be equal in both grids',...
                mfilename);
        end
    case 3
        if grOld.layersAreUniform
            if abs(grOld.zGr(  1) - grNew.zGr(  1))>delta || ...
               abs(grOld.zGr(end) - grNew.zGr(end))>delta
                error('%s: The bottom of the old grid (%g) and the new grid (%g) must match',...
                    mfilename,grOld.zGr(end),grNew.zGr(end));
            end
        else
            if any(abs(grOld.zGr(:,:,  1) - grNew.zGr(:,:,  1))>delta) || ...
               any(abs(grOld.zGr(:,:,end) - grNew.zGr(:,:,end))>delta)
                error('%s: The bottom of the old grid (%g) and the new grid (%g) must match',...
                    mfilename,grOld.zGr(end),grNew.zGr(end));
            end
        end
end

%% The to grid and the to values
if strcmpi(code,'zlist')
    % one value per layer, don't need to transfer entire grid, which may be large
    
    % elevation vector grid average
    
    wFr = grOld.zGr; wFr=unique(wFr(:)*ROUND)/ROUND;
    wTo = grNew.zGr; wTo=unique(wTo(:)*ROUND)/ROUND;
    
    wInt = unique([wFr; wTo]);
    wmInt=0.5*(wInt(1:end-1)+wInt(2:end));
    
    varTo = zeros(grNew.Nz,1);
    
    % indices from intermediate grid wInt into fr and to grid
    Ifr = min(numel(wFr)-1,floor(interp1(wFr,1:length(wFr),wmInt))); % index grid 1 in merged, size(Nint,1)
    Ito = min(numel(wTo)-1,floor(interp1(wTo,1:length(wTo),wmInt))); % index grid 2 in merged, size(Nint,1)
    
    % pull values from from grid into intermediategrid
    vInt = varFr(Ifr);
    
    % find start and end index of if to grid in intermediate grid
    [~,First] = unique(Ito,'first');  % First where they start in vInt
    [~,Last ] = unique(Ito,'last');   % and Last where they end
    
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
switch dim
    case {1,'y'}
        dFr  = reshape(permute(grOld.DY,[1,2,3]), [prod(grOld.size),1]);
        dTo  = reshape(permute(grNew.DY,[1,2,3]), [prod(grNew.size),1]);
        if ~strcmpi(code,'stress')
            varFr=     permute(varFr,   [1,2,3]);
        end
    case {2,'x'}
        dFr  = reshape(permute(grOld.DX,[2,1,3]), [prod(grOld.size),1]);
        dTo  = reshape(permute(grNew.DX,[2,1,3]), [prod(grNew.size),1]);
        if ~strcmpi(code,'stress')
            varFr=     permute(varFr,   [2,1,3]);
        end
    case {3,'z'}
        dFr  = reshape(permute(grOld.DZ,[3,2,1]), [prod(grOld.size),1]);
        dTo  = reshape(permute(grNew.DZ,[3,2,1]), [prod(grNew.size),1]);
        if ~strcmpi(code,'stress')
            varFr=     permute(varFr,   [3,2,1]);
        end
end

if any(dFr<grOld.MINDZ/2);
    error('%s: repeated elevations in dFr, layers in grOld must have thickness >= gr.MINDZ=%g\n',mfilename,grOld.MINDZ);
end
if any(dTo<grNew.MINDZ/2);
    error('%s: repeated elevations in dTo, layers in grNew must have thickness >= gr.MINDZ=%g\n',mfilename,grNew.MINDZ);
end

varFr=varFr(:);
varTo= zeros(size(dTo));

wFr  = unique(round([0; cumsum(dFr(:))]*ROUND)/ROUND); % from coordinates (L) in a single vector
wTo  = unique(round([0; cumsum(dTo(:))]*ROUND)/ROUND); % to   coordinates (L) in a single vector

%% Intermediate grid and intermediate values
wInt  = unique([wFr;wTo]);                 % intermediate grid-line coords
wmInt = 0.5*(wInt(1:end-1)+wInt(2:end));   % intermediate center coords
dInt =  diff(wInt);                        % intermediate grid cell widths
vInt  = zeros(size(dInt));                 % intermediate grid values

%% Pointers from intermediate grid into from grid Ifr and into target grid Ito
Ifr = min(numel(wFr)-1,floor(interp1(wFr,1:length(wFr),wmInt))); % index grid 1 in merged, size(Nint,1)
Ito = min(numel(wTo)-1,floor(interp1(wTo,1:length(wTo),wmInt))); % index grid 2 in merged, size(Nint,1)

%% Get unique varTo indices
[I,First] = unique(Ito,'first');  % First where they start in vInt
[~,Last ] = unique(Ito,'last');   % and Last where they end

%% Use only the overlapping grids

% this is implied by the coincidence of tops and bottoms of the to and the from grid. 
J  = find(~isnan(Ifr)); % J are valid indices in Ifr in case of non overlapping grid portions

if strcmpi(code,'stress')
    %% In this case the we deal with a stress list with LRC indices of which
    % the index for the current dimension has to be adapted to the new grid.
    % We pass the linear layer index. We rcompute it for the permuted
    % grid. The index then correspoinds with the
    % linear cell in the wFr vector. We can thus immediately compute
    % the center of the cell of the stress in the wFr vector. This
    % w-value (grid coordinate) is the same in the zTo vector, which is in
    % sequence of the idx of the new grid. We kan thus immediately look up the
    % linear index of the zNew array, compute the new layer and
    % transfer that back to the stress list.
    % The condition for this to work is that the top and
    % bottom planes of the from and the to networks are the same.

    % varFr contains the linear index of the cells
    switch dim
        case {1,'y'}
            CRL = cellIndices(varFr,grOld.size,'CRL');
            idx = cellIndex(CRL,[grOld.Ny,grOld.Nx,grOld.Nz]); 
            ymIdx = 0.5*(wFr(idx)+wFr(idx+1));
            varTo = min(floor(interp1(wTo,1:length(wTo),ymIdx)),length(wTo)); % linear index grid 2, size(Nint,1)
        case {2,'x'}
            RCL = cellIndices(varFr,grOld.size,'RCL');
            idx = cellIndex(RCL,[grOld.Ny,grOld.Nx,grOld.Nz]); 
            xmIdx = 0.5*(wFr(idx)+wFr(idx+1));
            varTo = min(floor(interp1(wTo,1:length(wTo),xmIdx)),length(wTo)); % linear index grid 2, size(Nint,1)
        case {3,'z'} % varFr = idxOld, varTo =idxNew
            LCR = cellIndices(varFr,grOld.size,'LCR');
            idx = cellIndex(LCR,[grOld.Nx,grOld.Nz,grOld.Ny]); 
            zmIdx = 0.5*(wFr(idx)+wFr(idx+1));
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
        varTo = permute(reshape(varTo,[grNew.Ny,grNew.Nx,grNew.Nz]),[1 2 3]);
    case 2
        varTo = permute(reshape(varTo,[grNew.Nx,grNew.Ny,grNew.Nz]),[2 1 3]);
    case 3
        varTo = permute(reshape(varTo,[grNew.Nz,grNew.Nx,grNew.Ny]),[3 2 1]);
end


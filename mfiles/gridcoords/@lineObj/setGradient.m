function O = setGradient(o,varargin)
% lineObj/setGradient --- to fix the head gradient at an outflowing boundary.
% The line object is adapted to one per layer which, together generated a
% fixed gradient outflowing boundary. This is done by placing ddrais at the
% bottom of the cells along this boundary in every layer and setting the
% specific conductance to gradient times k, where k is special, namely k =
% HK(iL)-HK(iL-1) for all layers except the last, where it isHK(iL).
%
% USAGE: newlineObj = lineObj.setGradient(gr,gradient,HK)
%
% newLineObj is size (lineObj(:),NLay) of lineObjects, each original
% lineObj is replaced by NLay lineObjects in the respective layers.
% The type of the original lineObj is immaterial, it will be changed to
% 'DRN' for the newLineObj. The old ones stay intact and should be deleted
% or replaced by the new ones to prevent multiple objects on the same
% boundary.
%
% gradient = desired head gradient across lineObj
%            a positive gradient is outward
%            for an inward given gradient use a negative value
% HK       = 3D array of conductivities (horizontal)
%
% ISSUES  notice that ther is not always a solution possible with given
% gradients. Also, the procedure requires that the LAYCON>0 so that the
% model has a free water table and the water table does not intersect with
% the top of the model. In that case, the gradient will be somewhat
% different from what was specified.
%
% TO 141105

    [gr,varargin] = getType(varargin,'gridObj',[]);
    if isempty(gr), error('%s requires gridObj as argument'); end
    [gradient,varargin] = getNext(varargin,'double',[]);
    if isempty(gradient) || ~isscalar(gradient)
        error('Second argument must be the gradient, i.e. numeric and scalar');
    end
    [HK, ~ ] = getNext(varargin,'double',[]);
    if isempty(HK) || ~all(size(HK)==gr.size)
        error('Third argument must be the HK array, i.e. numeric and size [%d,%d,%d]',gr.size);
    end
    
    O(numel(o),gr.Nlay) = lineObj();

    for io=numel(o):-1:1
        
        % Set type to 'DRN'
        [o(io).type] = deal('DRN');
        
        % Consider the vertical plane through the line object, which can be
        % regarded as a kind of "wall" for which we compute the correct
        % conductance per m length to be used so that the boundary will
        % imlement a fixed given hydraulic gradient. (If water table is
        % below the top of the aquifer system and the flow is outward.
        
        % get the Iy and Ix of the cells along the wall        
        CRL = cellIndices(o(io).Idx,gr.size,'CRL');
 
        % cell IDX for all layers at this wall        
        IDX = zeros(gr.Nz,numel(o(io).Idx));
        for iL=1:gr.Nlay
            CRL(:,end) = iL;
            IDX(iL,:) = cellIndex(CRL,gr);
        end
        
        %% Get conductivity along teh wall
        K = HK(IDX);                        % conductivities in wall
        
        % Then compute the conductivities required for the fixed gradient
        % conductances
        K = diff(K(end:-1:1,:));          % diff
        K = [K(end:-1:1,:) ; HK(IDX(end,:))]; % add original k at bottom
        
        % Compute the fixed-gradient conductanes (all specific per m
        % length, so not multiplied by the o(io).P.L.
        C = gradient * K;             % conductance that generates const.grad
        
        % Elevation of drains will be the bottom of the layer
        H = gr.ZBlay(IDX);                  % elevation is set to bottom of layers

        % replace vertexHdrs with z by dummies
        izCol = strmatchi('z',o(io).vertexHdr);
        for i=1:numel(izCol)
            o(io).vertexHdr{izCol} = ['obsolete' o(io).vertexHdr{izCol}];
        end
        ihCol = strmatchi('H',o(io).vertexHdr,'exact');
        if ihCol            
            for i=1:numel(ihCol)
                o(io).vertexHdr{ihCol} = ['obsolete' o(io).vertexHdr{ihCol}];
            end
        end
        
        % add new column with relative elevation
        
        o(io).vertex(      :,end+[1 2]) = 0;
        o(io).cellLineVals(:,end+[1 2]) = 0;

        vertHdr   = [o(io).vertexHdr {'zRel','H'}];
        
        izRel = strmatchi('zRel',vertHdr);
        iH    = strmatchi('H'   ,vertHdr,'exact');
        iC    = strmatchi('C'   ,vertHdr,'exact');

        o(io).vertexHdr = vertHdr;
        
        % Generate new lineObj per layer to replace the original one
        O(io,:) = o(io); % Copy old one to all layers (new ones)        

        sLine_ = o(io).sLine;
        sLine_([1 end]) = o(io).sCell([1 end]);

        for iLay=gr.Nlay:-1:1
            
            O(io,iLay).name = sprintf('%s_layer_%d',o(io).name,iLay);
            
            O(io,iLay).vertex(:,izRel) = iLay;
            O(io,iLay).vertex(:,iC   ) = interp1(o(io).sCell,C(iLay,:)',sLine_);
            O(io,iLay).vertex(:,iH   ) = interp1(o(io).sCell,H(iLay,:)',sLine_);

            O(io,iLay).cellLineVals(:,izRel) = iLay;
            O(io,iLay).cellLineVals(:,iC)    = C(iLay,:)'; 
            O(io,iLay).cellLineVals(:,iH)    = H(iLay,:)';
            O(io,iLay).Idx                   = IDX(iLay,:);
            
            for idx = 1:numel(o(io).Idx)
                O(io,iLay).P(idx).z   = gr.ZBlay(IDX(iLay,idx));
                O(io,iLay).P(idx).zR  = iLay;
                O(io,iLay).P(idx).iz  = iLay;
                O(io,iLay).P(idx).iLay= iLay;
                O(io,iLay).P(idx).idx = IDX(iLay,idx); 
            end
            O(io,iLay).V = {0 1};
        end
    end

end
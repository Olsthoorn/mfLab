function o = removeCBD(o)
%  Model = Model.removeCBD(gr) -- removes all CBD-layers (Confining BeD) from model
% thereby replacing VKCB of these layers by HK and VK respectively and by
% using the cell values for the new layers for all other parameters.
% o is an array of  model variables each stored in a ModelObj with fields
% name, var, type, where name is the name of the variable (HK, TRAN etc),
% var contains the variable itself, and type is its mfLab type (see
% Model.type;
%
% TO 120810

% first get the grid contained in the array of model variables (o) using
% method grid:
grOld = o.grid;

if all(grOld.LAYCBD == 0)
    return;
end

% Set LAYCBDnew equal to all zeros accoring to:
LAYCBDnew = zeros(length(grOld.LAYCBD)+sum(grOld.LAYCBD>0),1);

% generate the new grid
grNew = gridObj(grOld.xGr,grOld.yGr,grOld.Z,LAYCBDnew,grOld.MINDZ,grOld.AXIAL);

[L2N, C2N  LC2N, N2L] = convertLAYCBD(grOld.LAYCBD,LAYCBDnew);
% The natural layer sequence has total length Nlay+Nbd
% The input layer sequence is just 1:Nlay
% The input cbd   sequence is just 1:Ncbd
% With N natural layer stack, L the input layer stack and C the input CBD stack:
% L2N  input layer sequence to natural layer sequence:  N(.., L2N(:,2) =L(..,L2N(:,1)) = L;
% LC2N inp. lay. seq. to slots for cbd in nat. sequence N(..,LC2N(:,2) =L(..LC2N(:,1))
% CC2N input cbd   sequence to natural layer sequence   N(.., C2N(:,2) =C(..,C2N(:,1)) = C;

    for i =1:length(o)
        switch o(i).type
            case 'zlist'
                if strcmp(o(i).name,'LAYVKA'), LAYVKA=o(i).var; end % set apart need it later with VKA
                if strcmp(o(i).name,'VKA')   , VKA   =XS(o(i).var(:)); end % set apart as 3D, need it later with LAYVKA
                % put lay in natural stack
                Natural = NaN(grOld.Nlay+grOld.Ncbd,1);
                Natural(L2N( :,2),1) = o(i).var;
                % put corresponding lay value in cbd slots of natural stack
                Natural(C2N(:,2)) = o(i).var(C2N(:,1));
                o(i).var = Natural;
            case '3Dtime',
                % Do nothing; doesn't require any change
            case '3Dcbd'
                % Do nothing, see 3Dlay, needs HK, VK or TRAN, HY
            case '3Dlay'
                if strcmpi(o(i).name,'Z'), % skip Z, does not change when confining beds are removed
                    continue;
                else
                    if strcmp(o(i).name,'VKA'), VKA = o(i).var; end % need it later with LAYVKA
                    % Natural(:,:,natural)=Inputlayers)
                    Natural = NaN(grOld.Ny,grOld.Nx,grOld.Nlay+grOld.Ncbd);
                    Natural(:,:,L2N(:,2))=o(i).var;
                    Natural(:,:,LC2N(:,2))=Natural(:,:,LC2N(:,1));
                    o(i).var = Natural; % skip cbd, is done afterwards when all o have been processed
                end
            case 'struct'
                Natural = NaN(grOld.Ny,grOld.Nx,grOld.Nlay+grOld.Ncbd);
                if strmatchi('values',fieldnames(o(i).var));
                    for iper=1:length(o(i).var)
                        % fill layer slots with entire layer array
                        Natural(:,:,L2N( :,2))=o(i).var(iper).values;
                        % cbd slots in natural copy from lay
                        Natural(:,:,LC2N(:,2))=o(i).var(iper).values(:,:,LC2N(:,1));
                        % immediately replace old values array ith new one 
                        o(i).var = Natural;
                    end
                elseif strmatchi('term',fieldnames(o(i).var))
                    for iper=1:length(o(i).var);
                        for iterm=1:length(o(i).var(iper).term)
                            % fill layer slots with entire layer array
                            Natural(:,:,L2N( :,2))=o(i).var(iper).term{iterm}(:,:,L2N( :,1));
                            % all cbd slots in term are zero
                            Natural(:,:,LC2N(:,2))= 0;
                            % immediately replace old term with new one
                            o(i).var(iper).term{iterm} = Natural;
                            % cbd slots in Natural keep zero for budget terms
                        end
                    end
                else
                    error('%s: Struct with unknown field <<%s>>.',mfilename,o(i).type);
                end

            case 'stress'
                %% This works only if LAYCBDnew is all zeros (as is the case for this function)
                % replace the Z-indices by the new layer indices.
                o(i).var(:,2) = N2L(L2N(o(i).var(:,2),2),1);

            case {'wellObj','kDObj'}
                well = o(i).var;
                for iw = 1:length(well)
                    well(iw).iLay     = N2L(L2N(well(iw).iLay,2),1);
                    well(iw).LRC(:,1) = well(iw).iLay;
                    well(iw).idx      = cellIndex(well(iw).ix,well(iw).iy,well(iw).iLay,grNew.size);
                end
                o(i).var = well;
            case 'wellSeriesObj'
                % Do Nothing
            case 'gridObj'
                o(i).var = grNew;
            otherwise
                error('mfLab:removeCBD:unknownVariableType',...
                    '%s: Unknown type field <<%s>> for variable <<%s>>',mfilename,o(i).type,o(i).name);
        end
    end
    
    %% Layer Property Flow package uses VKA for vertical conductivity or vertical
    % anisotropy for the model layers depending on LAYVKA being zero or not.
    % It uses VKCB for the vertical conductivity of possible confining beds
    % (according to LAYCBD)
    iHK     = strmatchi('HK'    ,{o.name},'exact');
    iVKCB   = strmatchi('VKCB'  ,{o.name},'exact');
    iVKA    = strmatchi('VKA'   ,{o.name},'exact');
    iLAYVKA = strmatchi('LAYVKA',{o.name},'exact');
    iVK     = strmatchi('VK'    ,{o.name},'exact');
    
    if iVKCB % exists, put input CBD in CBD slots of natural stack in this case o(iHK).var
             % which as already been processed and has no CFB ans so
             % Nlay + Ncbd slots (layers)
        if iHK
            o(iHK).var(:,:, C2N(:,2)) = o(iVKCB).var;
        else     % put corresponding HK layers in CBD slots of natural stack 
            o(iHK).var(:,:,LC2N(:,2)) = o(iHK).var(:,:,LC2N(:,1));
        end
        if iVK
            o(iVK).var(:,:, C2N(:,2)) = o(iVKCB).var;
        else     % put corresponding HK layers in CBD slots of natural stack
           % o(iVK).var = o(iHK).var;
           % o(iVK).var(:,:,LC2N(:,2)) = o(iVK).var(:,:,LC2N(:,1));
        end
    end
    
    if iHK && ~iVK   % HK exists but VK not
        if ~(iVKA && iLAYVKA)
            error('%s: no variable VK is present so we need VKA, LAYVKA and HK, but not all of them found in Model',mfilename);
        else
            VK = o(iHK); % we reserved VKA above.
            VK.name = 'VK';
            
            for i = 1:LAYVKA
                if LAYVKA(i)>0,  % VKA is vertical anisotropy
                    VK.var(:,:,L2N(i,2)) = o(iHK).var(:,:,L2N(i,2))./VKA(:,:,i); % VKA stored before;
                else                     % VKA is vertical conductivity VK
                    VK.var(:,:,L2N(i,2)) = VKA(:,:,L2N(i,2));
                end
            end
            
            iVKCB   = strmatchi('VKCB',{o.name},'exact');
            if iVKCB % VKCB is given
                VK.var(:,:, C2N(:,2)) = o(iVKCB).var;
            else % VKCB is not given, use VK of overlying layer instead
                VK.var(:,:,LC2N(:,2)) = VK.var(:,:,LC2N(:,1));
            end
            o(iHK).var(:,:,C2N(:,2)) = VK.var(:,:,C2N(:,2));
        end
        o(end+1) = VK;
        if iVKCB,   o(iVKCB) = []; clear('iVKCB');   end
        if iLAYVKA, o(iLAYVKA)=[]; clear('iLAYVKA'); end
        if iVKA,    o(iVKA)   =[]; clear('iVKA');    end
    end

    % The Block-centered flow package uses TRAN for layer transmissivities
    % of layers that are not convertible (LAYCON 0 or 2) and VCONT adnd HY
    % for layers that are convertible (LAYCON 1 oe 3). It requires VCONT
    % specifying the vertical specific conductance between two model
    % layers.
    % mfLab aways requires the full TRAN array
    iHY   = strmatchi('HY'   ,{o.name},'exact');
    iTRAN = strmatchi('TRAN' ,{o.name},'exact');
    iVCONT= strmatchi('VCONT',{o.name},'exact');
    
    if iTRAN % obligatory for BCF package
        if iVCONT
            %Check formula: VCONT = vk/D -> vk = VCONT*D and TRAN = VCONT*D^2
            o(iTRAN).var(:,:,C2N(:,2)) = o(iVCONT).*gr.DZlay.*gr.DZlay;
        else
            % Check formula klay = Tlay/Dlay; Tcbd = Tlay/Dlay * Dcbd
            % Vertical anisotropy is implicit in VCONT and tranfers now to
            % the horizontal transmissivity of the 
            o(iTRAN).var(:,:,LC2N(:,2))= gr.DZcbd.*o(iTRAN).var(:,:,LC2N(:,1))./o(iTRAN).var(:,:,LC2N(:,1));
        end
    end
    if iHY
        if iVCONT
            % Check formula
            o(iHY).var(:,:,C2N(:,2)) = o(iVCONT).var./gr.DZcbd;
        else
            % Check formula
            o(iHY).var(:,:,LC2N(:,2)) = gr.DZcbd./o(iHY).var(:,:,LC2N(:,1));
        end
    end
    if iVCONT, o(iVCONT) = []; end
end


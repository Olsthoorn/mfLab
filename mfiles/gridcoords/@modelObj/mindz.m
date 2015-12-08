function o = mindz(o,MINDZ,layers)
% Model = Model.mindz(MINDZ,layers) -- sets minimum layer thickness while
% maintaining the transmissivity and vertical resistance of the
% combination. A mindz is important to reduce computation time in mt3dms
% and SEAWAT
%
% The minimum thickness is asserted at the bottom of the layer.
% 
% TO 120612

gr = o.grid();

% verify we have no confining beds
if any(gr.LAYCBD)~=0
    error(['%s: Sorry, byt this procdure is only valid if all LAYCBD==0,\n',...
          'first remove your CBD layers with gridOb.removeCBD'],mfilename);
end

% decide how MINDZ and layers are interpreted from the input
if nargin<3
    if size(MINDZ,1)==2,     % MINDZ =[layers ; mindz]
        layers = MINDZ(1,:);
        MINDZ  = MINDZ(2,:);
    else
        layers=1:gr.Nz;
        MINDZ=MINDZ(1)*ones(size(layers));
    end
else
    MINDZ=MINDZ(1)*ones(size(layers));
end

if any(MINDZ)<0, error('%s: MINDZ=%g must be greater than zero',mfilename); end

if any(layers<1) || any(layers>gr.Nz)
    error('%s: Layers must be from 1 to gr.Nz (%d)',gr.Nz);
end

% Get condutivities van Model array:
HK = o(strmatchi('HK',{o.name},'exact')).var;
VK = o(strmatchi('VK',{o.name},'exact')).var;

if ~all(size(HK)==gr.size), error(errstr,mfilename,'HK',size(HK),gr.size); end
if ~all(size(VK)==gr.size), error(errstr,mfilename,'VK',size(VK),gr.size); end

N  = gr.Nx*gr.Ny; % number of cells in layer to compute global indices
Z  = gr.Z;        % original elevations from grid

% itp (top plane == layer number)
[layers,I] = sort(layers,'ascend');
 MINDZ     = MINDZ(I);
 
 % layer number corresponds with that of its ceiling plane
for i=1:numel(layers)

    itp = layers(i);    

    % recursively from below as more than one layer may have to be lowered
    % for every layer iz we start checking from the bottom

    % bottom planes to check, ignoring the bottom plane of the model
    for ibp=gr.Nz:-1:itp+1

        % thickness between plane itp and ibp
        DZ = Z(:,:,itp)-Z(:,:,ibp);

        % where is it less than MINDZ?
        I = find(DZ<MINDZ(i));

        if isempty(I)
            % done, move check bottom plane one layer up
        else
            % Downward correction necessary for bottom plane
            dz   = MINDZ(i) - DZ(I);

            Itop = I+(ibp-2)*N; % global index of top of current layer where too thin 
            Ibot = I+(ibp-1)*N; % global index of bot of current layer where too thin

            if ibp<gr.Nz+1  % Don't shift the bottom of the model
                
                % Make sure that the new HK for this part of the layers is unaltered
                % Note that DZ(I)+dz equals MINDZ over I, length(dz)==length(I)
                % This sets the new HK(Itop) such that the tranmissivity over
                % DZ(I)+dz is unaltered
                Dtop      = max(gr.MINDZ,Z(Itop)-Z(Ibot));
                HK(Itop)  = (HK(Itop).* Dtop + HK(Ibot).*dz)./(Dtop+dz);

                % Make sure that the new VK for this part of the layer is unchanged
                % Note that (DZ(I)./ VK(Itop) + dz./VK(Itop+N) is total
                % vertical hydraulic resistance over DZ(I)+dz
                VK(Itop)  = (Dtop+dz)./(Dtop./ VK(Itop) + dz./VK(Ibot));

                % The correction of the underlying layer is implied by the change of
                % the vertical position of the bottom of the current layer,  i.e. the
                % top of the underlying layer
                Z (Ibot) = Z(Ibot) - dz;
            else % underlying layer is absent for last layer
                %% In fact, we don't lower the bottom of the model
                %HK(Itop)  = (HK(Itop).* Dthis)/(Dthis+dz);
                %VK(Itop)  = (Dthis+dz)./(Dthis./ VK(Itop));
                % Z (Ibot) = Z(Ibot) - dz; % exceeds gr.Nz+1
            end

        end
    end
end

%% correct zero thickness layers with gr.MINDZ
for itp=1:gr.Nz
    Z(:,:,itp+1) = min(Z(:,:,itp+1),Z(:,:,itp)-gr.MINDZ);
end

grNew = gridObj(gr.xGr,gr.yGr,Z,gr.LAYCBD,gr.MINDZ,gr.AXIAL);

o(strmatchi('gridObj',{o.type},'exact')).var = grNew;
o(strmatchi('HK'     ,{o.name},'exact')).var = HK;
o(strmatchi('VK'     ,{o.name},'exact')).var = VK;

%% Show old and new layer thickness:
% fprintf('%s:\n',mfilename);
% fprintf('Previous minimum layer thickness and new ones per layer, changed layers: ');
% fprintf('%5s %12s %12s\n','layer','start','end');
% fprintf('%5d %12g %12g\n',[(1:gr.Nz)',XS(min(min(gr.DZ,[],2),[],1)),XS(min(min(grNew.DZ,[],2),[],1))]');
% fprintf('\n');


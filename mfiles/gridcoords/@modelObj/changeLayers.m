function o = changeLayers(o,planes,subdivisions,dim)
% Model = Model.changeLayers(planes,subdivisions) -- replaces the layers of grid by new ones
% The new grid is completely resampled from the old using the new layers.
% o is a modelObj array, containing the arrays constituting the model
% gr is a gridObj
% planes      = layer planes to be kept
% subdivision = the subdivision of the vertical space between subsequent
%   new layers. Values are all integers >= 1.
% LAYCBD must be all zeros. Use modelObj.removeCBD() to achieve this if
% necessary.
%
% TO 120821

% First layer plane is implied, so remove if present in planes
if planes(1) == 1, planes = planes(2:end); end

% first get the grid contained in the array of model variables (o) using
% method grid:
grOld = o.grid;

if nargin<4, dim = 3; end

% set LAYCBDnew equal to all zeros accoring to:
if ~all(grOld.LAYCBD==0), error('%s: LAYCBD must be all zeros to use modelObj.changeLayers',mfilename); end

% generate the new grid
if isempty(planes), planes=(2:grOld.Nz+1); end

Znew = grOld.newZ(planes,subdivisions);

% new grid
grNew = gridObj(grOld.xGr,grOld.yGr,Znew,0,grOld.MINDZ,grOld.AXIAL);

%% make sure wellObj is the last item in o to deal with, because it needs HK
if exist('well','var')
    o = [o(~ismember({o.type},'wellObj')),o( ismember({o.type},'wellObj'))];
    o = [o(~ismember({o.type},'MNW1Obj')),o( ismember({o.type},'MNW1Obj'))];
    o = [o(~ismember({o.type},'MNW2Obj')),o( ismember({o.type},'MNW2Obj'))];
end
%% loop:
for i =1:length(o)
    fprintf('%s: converting variable %10s of type %10s\n',mfilename,o(i).name,o(i).type);
    if isempty(o(i).var), continue; end
    switch upper(o(i).type)
        case 'ZLIST'
            % let's do this for the mean Z. Apply to mean of Znew.
            o(i).var = grOld.gridTransfer(o(i).var,grNew,'zlist','z');
        case '3DTIME',
            % Do nothing; doesn't require any change
        case '3DCBD'
            % Do nothing, see 3Dlay, needs HK, VK or TRAN, HY
            %error('%s: There can be no "3Dcbd" type variable in this function',mfilename);
        case '3DLAY'
            switch o(i).name
                case 'VK', code='harmonic';
                case 'Z',
                    o(i).var = grNew.Z;
                    continue;
                otherwise
                 code = 'geometric';
            end
            if iscell(o(i).var), % e.g. STCONC
                for icell=1:numel(o(i).var)
                    o(i).var{icell} = grOld.gridTransfer(o(i).var{icell},grNew,code,dim);
                    if ismember(o(i).name(1),'ijklmnIJKLMN') % integers? -> round off
                        o(i).var{icell} = round(o(i).var{icell});
                    end
                end
            else
                o(i).var = grOld.gridTransfer(o(i).var,grNew,code,dim);
                if ismember(o(i).name(1),'ijklmnIJKLMN') % integers? -> round off
                    o(i).var = round(o(i).var);
                end
            end
            % keep HKnew in memory to use it with wellObj
            if strcmp(o(i).name,'HK'), HKnew = o(i).var; end
        case 'STRUCT'
            if strmatchi('values',fieldnames(o(i).var(1)));
                for iper=1:length(o(i).var)
                    o(i).var(iper).values = grOld.gridTransfer(o(i).var.values,grNew,'geometric',dim);
                    o(i).var(iper).NLay = size(o(i).var.values,3);
                end
            elseif strmatchi('term',fieldnames(o(i).var(1)))
                for iper=1:length(o(i).var);
                    for iterm=1:length(o(i).var(iper).term)
                        o(i).var(iper).term{iterm} = grOld.gridTransfer(o(i).var(iper).term{iterm},grNew,'geometric',dim);
                        o(i).var(iper).NLay = size(o(i).var(iper).term{iterm},3);
                    end
                end
            else
                error('%s: Struct with unknown field <<%s>>.',mfilename,o(i).type);
            end

        case 'STRESS'
            %% This works only if LAYCBDnew is all zeros (as is the case for this function)
            % replace the Z-indices by the new layer indices.
            % Procedure, put the stress in the center of the old cells
            % then interpolate the new cell numbers in the new grid based
            % on the elevations. Get the true layer numbers by subtracting
            % a multiple of the new number of layers
            if iscell(o(i).var)
                for ic=1:numel(o(i).var)
                    IdxOld =  cellIndex(o(i).var{ic}(:,[4,3,2]),grOld.size);
                    IdxNew =  grOld.gridTransfer(IdxOld,grNew,'stress',dim);
                    LRC = cellIndices(IdxNew,grNew.size,'LRC');
                    o(i).var{ic}(:,2) = LRC(:,1);
                end
            else
                IdxOld =  cellIndex(o(i).var(:,[4,3,2]),grOld.size);
                IdxNew =  grOld.gridTransfer(IdxOld,grNew,'stress',dim);
                LRC = cellIndices(IdxNew,grNew.size,'LRC');
                o(i).var(:,2) = LRC(:,1);
            end
        case {'WELLOBJ','KDOBJ','MNW1OBJ','MNW2OBJ'}
            % make sure wellObj comes last so that HKnew is present
            well = o(i).var;
            o(i).var = well.setWell(grNew,HKnew);
        case 'WELLSERIESOBJ'
            % Do Nothing
        case 'GRIDOBJ'
            o(i).var = grNew;
        otherwise
            error('mfLab:removeCBD:unknownVariableType',...
                '%s: Unknown type field <<%s>> for variable <<%s>>',mfilename,o(i).type,o(i).name);
    end
end


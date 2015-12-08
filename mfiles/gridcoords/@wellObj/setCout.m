function o=setCout(o,varargin)
%WELLOBJ/SETCOUT -- computes the output concentration of the wells
%
% USAGE:
%     well = well.setCout([iComp iComp2 ...][,B][,STCONC])
%     well = well.setCout(C1, C2, C3 ...[,B][,STCONC])
%     well = well.setCout(C ,iComp[,B][,STCONC])
%
% sets well(..).Cout, the computed well concentration. It is the flow-
% averaged concentration of the well's model cells during the simulation.
% It, therefore, serves as the well's output concentration. This is the flow
% averaged value. The flow is the total well flow times the fraction
% for each cell. Where the fraction is the transmissivity of the cell
% divided by the total transmissivity of the cells that take part in the well.
% If B, the struct read by readBud is specifed it will be used to mix
% the concentration at the well with the actual extractions, else
% well(iw).fQ will be use for tis purpose.
% Wells may also be multinode wells wells, for which the flow of each
% cell should be computed based on the computed cell-by-cell flows in B, not
% the cell-transmissivity derived fQ values.
% iComp is the compoment number in case the model compures more than one
% species. It will read files MT3D????.C, where ??? is a
% three didget iComp number. MT3D????.C files are produced by MT3DMS and
% SEAWAT.
% C, C1, C2 ... are structs read by readMT3D.m.
%
% if STCONC is given either as a numeric array of gr.size or a cell array
% where each element has size gr.size (for each species) the STCONC will be
% subtracted from the actual concentrations.
%
% Dissolved concentrations will be stored as well.Cout(iComp,1:NT)
%
% The well(iw).Q and injection conc well(iw).C that are necessary for MT3DMS
% and SEAWAT can be set using corresponding conlumns in the PER worksheet
% specified as arguments when calling the wellObj constructor; they don't
% need this function.
%
% SEE ALSO: wellObj
%
% TO 120423 121017 131125

[STCONC,varargin] = getType(varargin,'cell',{});
if isempty(STCONC)
    is = ~cellfun(@isscalar,varargin) & cellfun(@isnumeric,varargin);
    if any(is)
        STCONC = varargin(is);
        STCONC = STCONC(1);
        varargin(is)=[];
    end
end

% see if budget struct is provided
is = cellfun(@isstruct,varargin) & cellfun(@hasFieldTerm,varargin);
budget = any(is);
if budget
    varargin = [varargin(~is) varargin(is)];
end % put B at end
N = numel(find(~is));

% which class ??    
switch class(o)
    case 'wellObj', LBL = 'WELLS';
    case 'MNW1Obj', LBL = 'MNW'; % check label MODFLOW
    case 'MNW2Obj', LBL = 'MNW'; % check label MODFLOW
    otherwise
        error('%s: unknown class <<%s>>',mfilename,class(o));
end

% Just iComp(s) are specified in arbitrary order. We will read the UCN
% concentration files accordingly.

if all(cellfun(@isnumeric,varargin(1:N)))
    IComp= unique([varargin{1:N}],'descend');
    for iComp = IComp
            C{iComp}  = readMT3D(sprintf('MT3D%03d.UCN',iComp)); %#ok
    end
elseif all(cellfun(@isstruct,varargin(1:N)))
    IComp = 1:N;
    for iComp= N:-1:1
        C{iComp} = readMT3D(sprintf('MT3D%03d.UCN',iComp));
    end
elseif all(cellfun(@isstruct,varargin(1:2:N)) ...
       &   cellfun(@isscalar,varargin(2:2:N))) && rem(N,2)==0
    IComp = [varargin{2:2:end}];
    for i = numel(IComp):-1:1
        C{IComp(i)}=varargin{2*i-1};
    end
else
    error('unknown input argument combination');
end

%% Subtract STCONC if given
if ~isempty(STCONC)
    for iComp=IComp
        for it = numel(C{iComp}):-1:1
            C{iComp}(it).values = C{iComp}(it).values - STCONC{iComp};
        end
    end
end


% For all wells
for iw = numel(o):-1:1
    
    %% all components
    for iComp = numel(C):-1:1
        
        %% Dissolved species: if C for this component exists
        if ~isempty(C{iComp}) % if iComp was specified to be > 1
            
            %% All stress periods
            for it=numel(o(iw).Q):-1:1
                
                %% if wells are specified in this stress period
                if budget % we don't read in B for speed and memory, we use varargin{end} instead (call by ref)
                    iLbl = strmatchi(LBL,varargin{end}(it).label);
                    if iLbl
                        o(iw).Cout(iComp,it) = sum(C{iComp}(it).values(o(iw).idx(:)') .* ...
                            varargin{end}(it).term{iLbl}(o(iw).idx(:)')) / ...
                        sum(varargin{end}(it).term{iLbl}(o(iw).idx));
                    else
                        o(iw).Cout(iComp,it) = sum(C{iComp}(it).values(o(iw).idx(:)') .* ...
                            o(iw).fQ(:)');
                    end
                else
                    o(iw).Cout(iComp,it) = sum(C{iComp}(it).values(o(iw).idx(:)') .* ...
                        o(iw).fQ(:)');
                end
            end
        end        
    end
end
end

function hasit = hasFieldTerm(S)
  hasit = isfield(S,'term');
end

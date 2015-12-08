function [Model grNew] = mergeModel(Model,type,data,dim)
%MERGEMODEL convert Model into output Model with newCoords along dim
%
% Example:
%    [Model grNew] =gr.join(Model,newCoords,code)
%
%    real world coords are used for dim is 1 and 2 relative coords for dim 3
%    at least when z coordinates are nu uniform, in which case th newCoords
%    are relative z values starting with 0 at the top and declining downward
%
% These function are superseded by methods of the modelObj
%
% SEE ALSO: gridObj modelObj cutoutBCN cutoutXTA RefineBCN RefineGrid JoinBCN JoinLayers removeCBD
%
% TO 100601 120610

if ~ismember(dim,[1 2 3 'x' 'y' 'z'])
    error('%s: unknnown dimension should be one of [1 2 3 ''x'' ''y'' ''z'']',mfilename,dim);
end

if dim=='x',dim = 2; elseif dim=='y', dim=1; elseif dim=='z', dim=3; end

%% Merge the grid first and separately
iGRID = strmatchi('gridObj',{Model.type});

gr = Model(iGRID).var;

switch dim
    case 1, oldCoords    = gr.yGr(:)';
    case 2, oldCoords    = gr.xGr(:)';
    case 3, oldCoords    = gr.zGr(:)';
end

switch type
    case {'joinArray','mergeArray','merge','join'}, % intput is new relative coordinates
        newCoords = joinArray2Coords(data,oldCoords);
    case {'splitArray','refine','split'}, % input is a splitarray (suitable for layers)
        [newCoords, Ifirst, ~, Imid]= splitArray2Coords(data,oldCoords);
    case 'newcoords', newCoords = data;
end

grNew = mergeGrid(gr,newCoords,dim);

%% Merge the rest of of the data in Model
column = [3 4 2];  % the column(i) = dimension

for iM=1:numel(Model)
    switch Model(iM).type
        case 'zlist'
            if dim==3
                Model(iM).var=gridsTransfer(oldCoords,Model(iM).var,newCoords,'maximum',1);
            end
            
        case '3Dtime' % for instance rech or irich
            if ismember(dim,[1 2])
                if regexpi(Model(iM).name,'^[ijklmn]'),  code = 'abundant';
                else                                  code = 'geometric';
                end           
                Model(iM).var=gridsTransfer(oldCoords,Model(iM).var,newCoords,code,dim);
            end
        case 'stress'
            % We will only relink the cell to the joined cell or to one of
            % the new cells that are generated when cells are split. We
            % will not here generate new lines in the stress lists. So all
            % lines in the list remain and some may point to the same cell
            % after cells or layers have been joined.
            % It is probably the best compromise at this point.
            % The results depend on MODFLOW adding multiple stresses to the
            % same cell.
            switch type
                case 'join'
                    Model(iM).var(:,column(dim)) = data(2,Model(iM).var(:,column(dim)))';
                case 'split'
                    if dim == 3,
                        Model(iM).var(:,column(dim)) = Ifirst(Model(iM).var(:,column(dim)))';
                    else
                        Model(iM).var(:,column(dim)) = Imid(  Model(iM).var(:,column(dim)))';
                    end
            end
            
        case '3Dlay'
            if regexpi(Model(iM).name,'^v'),             code ='harmonic';
            elseif regexpi(Model(iM).name,'^[ijklmn]'),  code= 'median';
            else                                         code='geometric';
            end
            Model(iM).var=gridsTransfer(oldCoords,Model(iM).var,newCoords,code,dim);
            
        case 'struct'
            if strmatchi('values',fieldnames(Model(iM).var))
                for iper = 1:length(Model(iM).var)
                    Model(iM).var(iper).values=gridsTransfer(oldCoords,Model(iM).var(iper).values,newCoords,'geometric',dim);
                end
                
            elseif strmatchi('term',fieldnames(Model(iM).var))
                for iper=1:length(Model(iM).var)
                    for it = 1:length(Model(iM).var(iper).term)
                        Model(iM).var(iper).term{it}=gridsTransfer(oldCoords,Model(iM).var(iper).term{it},newCoords,'geometric',dim);
                    end
                end
                
            else
                error('%s: no fieldname values or term in Model{%d,2}',mfilename,iM);
            end
            
        case 'gridObj'
            Model(iGRID).var = grNew;
            
        case 'wellObj'
            fprintf('Processing wells ...'); tic;
            iHK = strmatchi('HK',{Model.name});
            if ~iHK,
                iHK = strmatchi('TRAN',{Model.name});
            end
            if ~iHK, error('%s: Model has no HK or TRAN',mfilename); end
            well=Model(iM).var;
            
            Model(iM).var = well.setWell(Model(iGRID).var,Model(iHK).var);
            
            toc; done;
            
        case 'wellSeriesObj'
            % Do nothing
            
        otherwise
            fprintf('%s: unknown variable type <<%s>>',mfilename,Model(iM).type);
    end

end

function uGrNew = joinArray2Coords(join,uGrOld)
    % Compute coordinates form joinArray
    % joinArray = [1 2 3 4 5 6 7 8 9 10 11 12 13 14; ...
    %              1 1 1 2 2 3 3 3 3  4  4  4  4  4]
    % is converted to relative coordinates:
    % uRelOld = [1 2 3 4 5 6 7 8 9 10 11 12 13 14]
    % uRelNew = [1     4   6       10          14]

    uRelOld = [join(1,:) join(1,end)+1];
    uRelNew = [join(1,1) find(diff(join(2,:),1,2))+1 join(1,end)+1];
    uGrNew    = interp1(uRelOld,uGrOld,uRelNew);
end

function [uGrNew, Ifirst, Ilast, Imid]= splitArray2Coords(split,uGrOld)
    %% Compute coordinates from split array
    % A split array only has the cells to be splitted
    % [ N1 N1 N3 N4;  n1 n2 n3 n4' m1 m2 m3 m4];

    %% split, divide and merge to get new coords
    N = numel(uGrOld)-1;

    %% Relative new coordinates implied by uGrOld an splitArray
    uRelNew{N} = 1; % store number of subcells in cell array
    for i=1:N, uRelNew{i}=1; end
    for i=1:size(split,2), uRelNew{split(1,i)}=ones(1,split(2,i))/split(2,i); end
    uRelNew = 1 + [0 cumsum([uRelNew{:}])];        
    uRelOld = 1:length(uGrOld);
    uGrNew  = interp1(uRelOld,uGrOld,uRelNew);
    uGrNew(end) = uGrOld(end);
    
    %% In case the cells are split, either new stress cells are required or
    %  a choice has to be made as to in which of the new cells the boundary 
    %  stress is going to be placed. In case dim==3, we choose the first
    %  i.e. the top cell. In other cases the center cel. With the Ilasat
    %  index vector also the last of each new cells can be selected to
    %  place the stress.
    
    %% Assemble splits:
    Isplit{N} = 1; for i=1:N, Isplit{i} = 1; end 
    for i=1:size(split,2), Isplit{split(1,i)}=zeros(1,split(2,i)); Isplit{split(1,i)}(1)=1; end
    Isplit = [ cumsum([Isplit{:}]); 1:N+sum(split(2,:)-1)];
    
    %% put stress left
    [~, Ifirst]= unique(Isplit(1,:),'first');
    [~, Ilast ]= unique(Isplit(1,:),'last');
    Imid = round(median([Ifirst; Ilast],1));
end

function gr = mergeGrid(gr,newCoords,dim)
    %   Anew = mergeGrid(gr,newCoords,dim)
    %   merge grid cells into direction dim according to new coords
    %
    % TO 120610
    
    permutes={[2 1 3];[1 2 3];[1 3 2]};
    
    switch dim
        % Note that gridsTransfer has x-orientation:
        case 1
            yGr = permute(newCoords,permutes{dim});
            
            Z  = gridsTransfer(gr.yGr,gr.Z,newCoords,'geometric',dim);
            
            gr = gridObj(gr.xGr,yGr,Z,gr.LAYCBD,gr.MINDZ,gr.AXIAL);
            
        case 2
            % Note that gridsTransfer works in x-direction, no permute
            % necessary here:
            xGr        = newCoords;
            
            Z = gridsTransfer(gr.xGr,gr.Z,newCoords,'geometric',dim);

            gr = gridObj(xGr,gr.yGr,Z,gr.LAYCBD,gr.MINDZ,gr.AXIAL);

        case {3,'z'}
            % Note that gridsTransfer works in x-direction:
            % In z-direction we only use relative coordinates. For this the
            % newCoords and the oldCoords are based on gr.zGr which is a
            % z-vector of the layer-average z-values.
            oldCoords = gr.zGr(:)';
            %% Option 'width' is only necessary for DZ
            DZ           = gridsTransfer(oldCoords,gr.DZ,newCoords,'width',dim);
            Z            = repmat(gr.Z(:,:,1),[1,1,length(newCoords)]);
            Z(:,:,2:end) = Z(:,:,2:end) - cumsum(DZ,3);            
            LAYCBD = 0; % gridObj will extend it
            
            gr = gridObj(gr.xGr,gr.yGr,Z,LAYCBD,gr.MINDZ,gr.AXIAL);
    end

end

end
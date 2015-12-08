function o = merge(o,type,data,dim)
%MERGE convert Model into output o with newCoords along dim
% is a method of modelObj
%
% Example:
%    Model = Model.merge(o,type,data,dim)
%
% Inputs
%    type is oneof { 'joinArray','mergeArray','merge','join',...
%                    'splitarray','refine','split'}
%    data is the join or split array
%    dim  is the dimension along which the merging or splitting is done (use
%       1,2,3 or 'y','x', or 'z')
%
% Real world coords are used for dim is 1 and 2 relative coords for dim 3
% at least when z coordinates are nu uniform, in which case th newCoords
% are relatie z values starting with 0 at the top and declining downward
%
% TO 100601 120610 120814

if ~ismember(dim,[1 2 3 'y' 'x' 'z'])
    error('%s: unknnown dimension should be one of [1 2 3 ''y'' ''x'' ''z'']',mfilename,dim);
end

% use Matlab dim (y first)
if dim=='x',dim = 2; elseif dim=='y', dim=1; elseif dim=='z', dim=3; end

%% Merge the grid first and separately

grOld = o.grid();

switch type
    case {'joinArray','mergeArray','merge','join'}, % intput is new relative coordinates
        grNew = joinArray2Coords(grOld,data,dim);
        
    case {'splitArray','refine','split'}, % input is a splitarray (suitable for layers)
       [grNew, Ifirst, ~, Imid]= splitArray2Coords(grOld,data,dim);
end

%% Merge the rest of of the data in o
column = [3 4 2];  % the column(i) = dimension in stress lists

for iM=1:numel(o)
    switch o(iM).type
        case 'zlist'
            if dim==3
                o(iM).var=grOld.gridTransfer(o(iM).var,grNew,'zlist',1);
            end
            
        case '3Dtime' % for instance rech or irich
            if dim==1 || dim ==2
                if regexpi(o(iM).name,'^[ijklmn]'),
                    code = 'abundant';
                else
                    code = 'geometric';
                end
                o(iM).var=gridTransfer(...
                    gridObj(grOld.xGr,grOld.yGr,-(0:size(o(iM).var,3)),grOld.LAYCBD,grOld.MINDZ,grOld.AXIAL),...
                    o(iM).var,...
                    gridObj(grNew.xGr,grNew.yGr,-(0:size(o(iM).var,3)),grNew.LAYCBD,grNew.MINDZ,grNew.AXIAL),...
                    code,dim);
            else
                % ignore
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
                    if iscell(o(iM).var)
                        for ic=1:numel(o(iM).var)
                            o(iM).var{ic}(:,column(dim)) = data(2,o(iM).var{ic}(:,column(dim)))';
                        end
                    else
                        o(iM).var(:,column(dim)) = data(2,o(iM).var(:,column(dim)))';
                    end
                case 'split'
                    if dim == 3, % z
                        if iscell(o(iM).var)
                            for ic=1:numel(o(iM).var)
                                o(iM).var{ic}(:,column(dim)) = Ifirst(o(iM).var{ic}(:,column(dim)))';
                            end
                        else
                            o(iM).var(:,column(dim)) = Ifirst(o(iM).var(:,column(dim)))';
                        end
                    else
                        if iscell(o(iM).var)
                            for ic=1:numel(o(iM).var)
                                o(iM).var{ic}(:,column(dim)) = Imid(  o(iM).var{ic}(:,column(dim)))';
                            end
                        else
                            o(iM).var(:,column(dim)) = Imid(  o(iM).var(:,column(dim)))';
                        end
                    end
            end
            
        case '3Dlay'
            %display(o(iM));
            %display(o(iM).var);
            if regexpi(o(iM).name,'^v'),             code = 'harmonic';
            elseif regexpi(o(iM).name,'^[ijklmn]'),  code = 'median';
            else                                     code = 'geometric';
            end
            o(iM).var = grOld.gridTransfer(o(iM).var,grNew,code,dim);
            %display(o(iM).var);
            %fprintf \n
            
        case 'struct'
            if strmatchi('values',fieldnames(o(iM).var))
                for iper = 1:length(o(iM).var)
                    o(iM).var(iper).values = grOld.gridTransfer(o(iM).var(iper).values,grNew,'geometric',dim);
                end
                
            elseif strmatchi('term',fieldnames(o(iM).var))
                for iper=1:length(o(iM).var)
                    for it = 1:length(o(iM).var(iper).term)
                        o(iM).var(iper).term{it} = grOld.gridTransfer(o(iM).var(iper).term{it},grNew,'geometric',dim);
                    end
                end
                
            else
                error('%s: no fieldname values or term in o(%d).var}',mfilename,iM);
            end
            
        case 'gridObj'
            
            o(strmatchi('gridObj',{o.type})).var = grNew;
            
        case 'wellObj'
            fprintf('Processing wells ...'); tic;
            iHK = strmatchi('HK',{o.name});
            if ~iHK,
                iHK = strmatchi('TRAN',{o.name});
            end
            if ~iHK, error('%s: o has no HK or TRAN',mfilename); end
            well=o(iM).var;
            
            o(iM).var = well.setWell(grNew,o(iHK).var);
            
            toc; done;
            
        case 'wellSeriesObj'
            % Do nothing
            
        otherwise
            error('%s: unknown variable type <<%s>>',mfilename,o(iM).type);
    end

end

fprintf \n

function grNew = joinArray2Coords(grOld,join,dim)
    % Compute coordinates from joinArray
    % joinArray = [1 2 3 4 5 6 7 8 9 10 11 12 13 14; ...
    %              1 1 1 2 2 3 3 3 3  4  4  4  4  4]
    % is converted to relative coordinates:
    % uRelOld = [1 2 3 4 5 6 7 8 9 10 11 12 13 14]
    % iGr = [1     4   6       10          14]

    if any(diff(join(1,:))~=1)
        error('%s:joinArray line1 must contain all layer numbers in sequence.',mfilename);
    end
    
    iGr = [join(1,1) find(diff(join(2,:),1,2))+1 join(1,end)+1];
    Inew=unique(join(2,:));
    switch dim
        case {1,'y'}
            Z = grOld.Z(Inew,:,:);
            yGr    = grOld.yGr(iGr);
            for iy= Inew
                Z(iy,:,:) = mean(grOld.Z(join(1,join(2,:)==iy),:,:),dim);
            end
            grNew  = gridObj(grOld.xGr,yGr,Z,grOld.LAYCBD,grOld.MINDZ,grOld.AXIAL);
        case {2,'x'}
            Z = grOld.Z(:,Inew,:);
            xGr    = grOld.xGr(iGr);
            for ix= Inew
                Z(:,ix,:) = mean(grOld.Z(:,join(1,join(2,:)==ix),:),dim);
            end
            grNew  = gridObj(xGr,grOld.yGr,Z,grOld.LAYCBD,grOld.MINDZ,grOld.AXIAL);
        case {3,'z'}
            grNew  = gridObj(grOld.xGr,grOld.yGr,grOld.Z(:,:,iGr),grOld.LAYCBD,grOld.MINDZ,grOld.AXIAL);
    end

end

function [grNew, Ifirst, Ilast, Imid]= splitArray2Coords(grOld,split,dim)
    %% Compute coordinates from split array
    % A split array only has the cells to be splitted
    % [ N1 N1 N3 N4;  n1 n2 n3 n4' m1 m2 m3 m4];

    Nold = grOld.size(dim);
    
    %% complete the splitarray
    SPLIT = ones(2,Nold);
    SPLIT(1,:)=1:Nold;
    SPLIT(2,split(1,:))=split(2,:);
    
    Nnew = sum(SPLIT(2,:));
    
    Ifr   = NaN(Nnew,1)  ; Ifr(1) = 1; % Ifrom -> Ifr(size(1,Nnew)

    %% Relative new coordinates implied by uGrOld an splitArray
    
    % store number of subcells in cell array
    % default is all 1
    switch dim
        case 1
            yGr = NaN(Nnew+1,1); yGr(1)=grOld.yGr(1);
            k=1;
            for i=1:Nold
                for j=1:SPLIT(2,i)
                    yGr(k+1) = yGr(k) - grOld.dy(i)/SPLIT(2,i);
                    Ifr(k)=i;
                    k=k+1;
                end
            end            
            grNew = gridObj(grOld.xGr, yGr,grOld.Z(Ifr,:,:),grOld.LAYCBD,grOld.MINDZ,grOld.AXIAL);
        case 2
            xGr = NaN(Nnew+1,1); xGr(1)=grOld.xGr(1);
            k=1;
            for i=1:Nold
                for j=1:SPLIT(2,i)
                    xGr(k+1) = xGr(k) + grOld.dx(i)/SPLIT(2,i);
                    Ifr(k)=i;
                    k=k+1;
                end
            end            
            grNew = gridObj(   xGr,grOld.yGr,grOld.Z(:,Ifr,:),grOld.LAYCBD,grOld.MINDZ,grOld.AXIAL);
        case 3
            Znew = NaN(grOld.Ny,grOld.Nx,Nnew+1); Znew(:,:,1)=grOld.Z(:,:,1);
            k=1;
            for i=1:Nold
                for j=1:SPLIT(2,i)
                    Znew(:,:,k+1) = Znew(:,:,k) - grOld.DZ(:,:,i)/SPLIT(2,i);
                    Ifr(k)=i;
                    k=k+1;
                end
            end            
            grNew = gridObj(grOld.xGr,grOld.yGr,Znew,grOld.LAYCBD,grOld.MINDZ,grOld.AXIAL);
    end
    
    
    %% In case the cells are split, either new stress cells are required or
    %  a choice has to be made as to in which of the new cells the boundary 
    %  stress is going to be placed. In case dim==3, we choose the first
    %  i.e. the top cell. In other cases the center cel. With the last
    %  index vector also the last of new cells can be selected to place it in the new stress.
    
    %% put stress left
    [~, Ifirst] = unique(Ifr,'first');
    [~, Ilast ] = unique(Ifr,'last');
    Imid = round(median([Ifirst Ilast],2));
end

end
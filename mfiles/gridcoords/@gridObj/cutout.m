function [model] = cutout(o,model,Ix,Iy,Iz)
%% [model] = gr.cutout(model,Ix,Iy,Iz)
%  cuts variable model.var of mfLab type model.type from the grid (o) within limiting indices Ix,Iy,Iz
%  called from Model.submodel() see ModelObj.
%  model is of class modelObj (a struct containing name var and type of a
%     model variable like HK, STRTHD, well, et c
%  The size of the new grid is [1:length(Iy),1:length(Ix,gr.Nz,1:length(Iz)]
%  Ix, Iy and Iz must be contiguous valid indices in the existing (old) grid
%
%  The integrety of Ix, Iy and Iz has been checked before the call, but we
%  check it here again to prevent bypassing it.
%
% TO 120606

if ~all(Ix >= 1 & Ix <= o.Nx), error('%s: Ix outside limits of source grid)'); end
if ~all(Iy >= 1 & Iy <= o.Ny), error('%s: Iy outside limits of source grid)'); end
if ~all(Iz >= 1 & Iz <= o.Nz), error('%s: Iz outside limits of source grid)'); end

Nx     = length(Ix);
Ny     = length(Iy);
Nlay   = length(Iz);

switch model.type
    case 'gridObj',
        LAYCBD = o.LAYCBD;
        iznew = reshape(cumsum(reshape([ones(size(LAYCBD)) LAYCBD]',[2*length(LAYCBD),1])),[2,length(LAYCBD)])';
        iznew(LAYCBD==0,2)=0;
        iznew = iznew(Iz,:);
        iznew = reshape(iznew',[2*Nlay,1]);
        iznew = iznew(iznew~=0);
        iznew = [iznew; iznew(end)+1];

        model.var = gridObj(o.xGr([Ix Ix(end)+1]),o.yGr([Iy; Iy(end)+1]),o.Z(Iy,Ix,iznew),...
              LAYCBD,o.MINDZ,o.AXIAL);
        % do nothing
    case 'wellSeriesObj'
        % do nothing
        
    case 'zlist'
        model.var = model.var(Iz);
        
    case {'wellObj','kDObj'}
        well=model.var;
        wellIx = NaN(numel(well,1));
        wellIy = NaN(numel(well,1));
        for iw = numel(well):-1:1
            wellIx(iw) = well(iw).ix(1);
            wellIy(iw) = well(iw).iy(1);
        end
        
        % select well within selected region
        well= well(wellIx>=Ix(1) & wellIx<=Ix(end) & wellIy>=Iy(1) & wellIy<=Iy(end));

        idxOk = false(1,length(well));
        for iw=1:length(well)
            well(iw).ix   = well(iw).ix   - Ix(1) + 1;
            well(iw).iy   = well(iw).iy   - Iy(1) + 1;
            
            well(iw).iLay = well(iw).iLay(well(iw).iLay>=Iz(1) & well(iw).iLay<=Iz(end))-Iz(1)+1;
            
            if ~isempty(well(iw).iLay)
                well(iw).idx = cellIndex(  well(iw).ix, well(iw).iy, well(iw).iLay,[Ny,Nx,Nlay]);
                well(iw).LRC = cellIndices(well(iw).idx,[Ny,Nx,Nlay],'LRC');
            else
                well(iw).idx=NaN;
                well(iw).LRC=NaN;
            end
            idxOk(iw) = ~isnan(well(iw).idx(1));
        end
        
        well = well(idxOk);
        
        model.var=well;

    case 'stress'
        model.var = model.var( model.var(:,2)>=min(Iz) & model.var(:,2)<=max(Iz),:);
        model.var = model.var( model.var(:,3)>=min(Iy) & model.var(:,3)<=max(Iy),:);
        model.var = model.var( model.var(:,4)>=min(Ix) & model.var(:,4)<=max(Ix),:);

        model.var(:,2) = model.var(:,2)-Iz(1)+1;
        model.var(:,3) = model.var(:,3)-Iy(1)+1;
        model.var(:,4) = model.var(:,4)-Ix(1)+1;

    case '3Dlay',
        if strcmp(model.name,'C'), Iz = Iz(Iz<=o.Nlay-1); end
        if iscell(model.var) % e.g. STCONC with more components
            for icell = 1:numel(model.var)
                model.var{icell}=model.var{icell}(Iy,Ix,Iz);
            end
        else
            model.var=model.var(Iy,Ix,Iz);
        end
        
    case '3Darray',
        model.var=model.var(Iy,Ix,Iz);
        
    case '3Dcbd',
        iNew = cumsum(o.LAYCBD);
        iNew(o.LAYCBD==0)=0;
        iNew = iNew(Iz);
        iNew = iNew(iNew~=0);
        model.var = model.var(Iy,Ix,iNew);
        warning('gridOjb:cutout:VARtypeNotImplemented',...
            '%s: cutout not implemented for vartype <<%s>>',mfilename,model.type);
        
    case '3Dtime'
        model.var = model.var(Iy,Ix,:);
        
    case 'struct'
        fldnms = fieldnames(model.var);
        for iper=1:length(model.var)
            for ifld =1:length(fldnms)
                switch fldnms{ifld}
                    case 'values',
                        model.var(iper).(fldnms{ifld}) = model.var(iper).(fldnms{ifld})(Iy,Ix,Iz);
                    case 'term'                    
                        for j=1:length(model.var(iper).(fldnms{ifld}))
                            model.var(iper).(fldnms{ifld}){j} = model.var(iper).(fldnms{ifld}){j}(Iy,Ix,Iz);
                        end
                    case 'NROW'
                        model.var(iper).NROW = length(Iy);
                    case 'NCOL'
                        model.var(iper).NCOL = length(Ix);
                    case 'NLAY'
                        model.var(iper).NLAY = length(Iz);
                    case 'cols'
                        model.var(iper).cols = 1:length(Ix);
                    case 'rows'
                        model.var(iper).rows = 1:length(Iy);
                    case 'lays'
                        model.var(iper).lays = 1:length(Iz);
                    otherwise
                        % Do nothing
                end
            end
        end
        
    otherwise
        error('gridObj:cutout:unknownVarType', ...
        ' %s: vartype <<%s>> unknown see help for legal types',mfilename,model.type);
end

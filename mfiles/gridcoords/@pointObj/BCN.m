function cellList = BCN(o,varargin)
    %BCN: is a method of pointObj, lineObj and area2Obj
    %
    % USAGE: cellList = obj.BCN(gr,type,basename)
    %
    % obj/BCN generates a list of stress boundary tuples to write
    % to the respective MODLFOW stress file for this stress type.
    % BCN stands for Boundary Condition or STRESS (WEL, DRN, GHB,
    % RIV etc) and MODFLOW requires per stress period a list of the
    % concerning cells with per line [Layer Row Col data ...]
    % where data depends on the stress type. See Modflow manual.
    % In preparation of that, mfLab in this function generates a
    % cell array with one list of such tuples for each stress
    % period. mfLab addes the stress period number as first column,
    % so that for each tuple the stress period is known to which it
    % belongs and it does not matter in which order mfLab generates
    % them. All lists with tuples from whatever object it is generated
    % will be unpacked and sorted before writing the MODFLOW stress file.
    % You can always show what finally entered the modflow stress
    % file by plotting the cells in it using the gridObj method
    % gridObj.showWhatIsInMyStressFile(ModflowStressFile). See the
    % help of that function for details.
    % Notice that parts of the object what are outside the model
    % grid (either above, below left right, frong or behind) are
    % kicked out. This can be done, because the grid information is
    % known and contained in the gridObj and the coordinates of the
    % vertices of the pointObj,lineObj or area2Obj are contained in
    % these objects.
    %
    % This function is largely independent of the subclases lineObj
    % and area2Obj.
    %
    % USAGE:  obj = obj.BCN(type,gr)
    %  where type can be 'NIL','FLUX','WEL','DRN','GHB','RIV','CHD'
    %  and gr is the gridObj and further obj is of class poitObj,
    %  lineObj or area2Obj.
    %
    % The grid is an obligatory argument and can be placed anywhere
    % in the argument list.
    %
    % SEE ALSO: gridObj, pointObj, lineObj area2Obj gridObj/showWhatIsInMyStressFile
    %
    % TO 131001

    [gr        ,varargin]  = getType(varargin,'gridObj',[]);
    [stressType,varargin]  = getNext(varargin,'char','');
    [basename,    ~     ]  = getNext(varargin,'char','');

    [perNams,perVals,Nper] = getPeriods(basename,'PER');

    % Make sure headers are unique
    NU = notUnique(perNams);
    if ~isempty(NU)
        error('Headers in PER sheet are non unique <<%s>>',...
            sprintf(' %s',NU{:})); %#ok
    end

    % Return silently if o is empty
    if isempty(o)
        cellList{Nper,1} = [];
        return
    end

    % Assert stressType is one of the legal stresses
    if ~strmatchi(stressType,o(1).stressNames)
        cellList{Nper,1} = [];
        return;
    end

    % Assert that we have at least one object of the specified
    % stress stressType, else, return silently
    if strcmpi('WEL',stressType)
        o = o(strmatchi({stressType,'FLUX'},{o.type}));
    else
        o = o(strmatchi(stressType,{o.type}));
    end            
    if  isempty(o)
        cellList{Nper,1} = [];
        return;
    end

    % Initialize cell list
    lastCell = 0; 
    % We will generate a list of matlab cells (a Matlab cell array),
    % with one cell for each active stress period
    % containing the modflow cells(one line per cell i..e
    %     { ip L R C values;
    %       ip L R C values;
    %       ...
    %     }
    %     { ip+1 L R C values;
    %       ip+1 L R C values;
    %       ....
    %     }
    %    etc.
    % pertaiing active in this stress period. (ip is the stress
    % period number).
    % The cell array will be used in writeBCN to print out the
    % Modflow stess file.

    msgId = 'BCN:belowObj';
    warning('on',msgId);

    % Process objects in turn as each one may be active in arbitrry
    % stress periods, not necessarily all.
    for io=numel(o):-1:1

        jx = o(io).iColx;
        jy = o(io).iColy;

        % See if a PER-sheet column was specified to indicate when the
        % objects are active.
        % If so, then assert immediately that that column exists in the PER sheet.
        if ~isempty(o(io).activePerColName)
            try
                if ismember(o(io).activePerColName(1),'~-')  % use opposite values
                    ACTIVE = ~perVals(:,strmatchi(o(io).activePerColName(2:end),perNams,'exact'));
                else
                    ACTIVE =  perVals(:,strmatchi(o(io).activePerColName(1:end),perNams,'exact'));
                end
                ACTIVE = ACTIVE ~= 0;
            catch ME
                if ~iActivePerCol
                    error(ME.messageId,...
                        'can''t find column <<%s>> in PER sheet to see when %s <<%s>> is active',...
                        o(io).activePerColName,o(io).class);
                end
            end
        else
            ACTIVE = (1:Nper)';
        end

        % Which columns to write to the MODFLOW input files ?
        % There area at most 3 columns to be written (RIV has
        % three, GHB, DRN, CHD have 2 and WEL|FLUX has 1)

        ICol = zeros(1,3) ; % column indices into vertexHdr, cellVals               
        IColT= zeros(1,3);  % column indices into vertexTxtHdr, holding
        Irel = false(1,3);  % interpret head and elevations as relative or absolute values

        % ICol  pertains to numeric constant with time data
        % IColT pertains to varying with time data. I.e. through a string that is interpreted
        % as the header of a column with transient data in the PER
        % worksheet.
        % If columns containing elevation data have a header
        % starting with H they are considered absolute elevations,
        % if their header starts with D they are considered
        % relative with respect to ground surface. There elevation
        % is groundSurface + their value.
        % We have already checked at in the object constructor that
        % the necessary columns are all present. We will omit this
        % checking here.

        switch stressType
            case {'WEL','FLUX'}
                ICol (1) = strmatchi({'WEL','FLUX','Q'},o(io).vertexHdr,'exact','once');
                if ~ICol(1)
                    IColT(1) = strmatchi({'WEL','FLUX','Q'},o(io).vertexTxtHdr,'exact','once');
                end
                if ~any([ICol(1) IColT(1)])
                    error('No FLUX|WEL column specified for %s',o(io).type);
                end
            case 'DRN'

                msgId = 'mfileName:DRNpointsAboveGround';
                warning('on',msgId);

                ICol (1) = strmatchi({'H','D'},o(io).vertexHdr,'exact','once');
                if ICol( 1)
                    Irel(1)  = ismember(o(io).vertexHdr{ICol(1)}(1),'D');
                else
                    IColT(1) = strmatchi({'H','D'},o(io).vertexTxtHdr,'exact','once');
                    Irel( 1) = ismember(o(io).vertexTxtHdr{IColT(1)}(1),'D');
                end

                ICol(     2) = strmatchi({'C'},o(io).vertexHdr,'exact','once');
                if ~ICol( 2)
                    IColT(2) = strmatchi({'C'},o(io).vertexTxtHdr,'exact','once');
                end

            case 'GHB'
                ICol (1) = strmatchi({'H','D'},o(io).vertexHdr,'exact','once');
                if ICol( 1)
                    Irel(1) = ismember(o(io).vertexHdr{ICol(1)}(1),'D');
                else
                    IColT(1) = strmatchi({'H','D'},o(io).vertexTxtHdr,'exact','once');
                    Irel( 1) = ismember(o(io).vertexTxtHdr{IColT(1)}(1),'D');
                end

                ICol (2) = strmatchi({'C'},o(io).vertexHdr);
                if ~ICol( 2)
                    IColT(2) = strmatchi({'C'},o(io).vertexTxtHdr);
                end

                if ~any([ICol(1) IColT(1)]) && ~any([ICol(2) IcolT(2)])
                    error('No head|stage column or C column specified for %s',o(io).type);
                end
            case 'RIV'
                ICol (1) = strmatchi({'H','D'},o(io).vertexHdr,'exact','once');
                if ICol( 1)
                    Irel( 1) = ismember(o(io).vertexHdr{ICol(1)}(1),'D');
                else
                    IColT(1) = strmatchi({'H','D'},o(io).vertexTxtHdr,'exact','once');
                    Irel( 1) = ismember(o(io).vertexTxtHdr{IColT(1)}(1),'D');
                end

                ICol (2) = strmatchi({'C'},o(io).vertexHdr,'exact');
                if ~ICol( 2)
                    IColT(2) = strmatchi({'C'},o(io).vertexTxtHdr,'exact');
                end

                ICol (3) = strmatchi({'Hb','Db'},o(io).vertexHdr,'once');
                if ICol( 3)
                    Irel( 3) = ismember(o(io).vertexHdr{ICol(3)}(1),'D');
                else
                    IColT(3) = strmatchi({'Hb','Db'},o(io).vertexTxtHdr,'once');
                    Irel( 3) = ismember(o(io).vertexTxtHdr{IColT(3)}(1),'D');
                end

                if ~any([ICol(1) IColT(1)]) && ~any([ICol(2) IColT(2)]) && ~any([ICol(3) IColT(3)])
                    error('No head column or C column or Hb oclumn specified for %s',o(io).type);
                end

            case 'CHD'
                ICol(1)  = strmatchi({'H1','Hst','D1','Dst'},o(io).vertexHdr,'exact','once');
                if ~ICol( 1)
                    IColT(1) = strmatchi({'H1','Hst','D1','Dst'},o(io).vertexTxtHdr,'exact','once');
                end
                if ~(ICol(1) || IColT(1))
                     ICol(1) = strmatchi({'H','D'},o(io).vertexHdr,'exact','once');
                    if ~ICol( 1)
                        IColT(1) = strmatchi({'H','D'},o(io).vertexTxtHdr,'exact','once');
                    end
                end
                if ICol( 1)
                    Irel(1) = ismember(o(io).vertexHdr{ICol(1)}(1),'D');
                else
                    Irel(1) = ismember(o(io).vertexTxtHdr{IColT(1)}(1),'D');
                end

                ICol(2)  = strmatchi({'H1','Hst','D1','Dst'},o(io).vertexHdr,'exact','once');
                if ~ICol( 2)
                    IColT(2) = strmatchi({'H1','Hst','D1','Dst'},o(io).vertexTxtHdr,'exact','once');
                end
                if ~(ICol(2) || IColT(2))
                     ICol(2) = strmatchi({'H','D'},o(io).vertexHdr,'exact','once');
                    if ~ICol( 2)
                        IColT(2) = strmatchi({'H','D'},o(io).vertexTxtHdr,'exact','once');
                    end
                end
                if ICol( 2)
                    Irel(2) = ismember(o(io).vertexHdr{ICol(2)}(1),'D');
                else
                    Irel(2) = ismember(o(io).vertexTxtHdr{IColT(2)}(1),'D');
                end                        

            otherwise
                error('unrecognized stress type <<%s>>',stressType);
        end

        % In case the values are a string, it refers to a
        % column of the same hdr in the PER sheet, get that column
        % number from the PER sheet:
        IColPer = zeros(size(IColT));

        for ic = 1:numel(IColT)
            if IColT(ic)
                str         = o(io).vertexTxt{1,IColT(ic)}; % first point of this column
                ii = strmatchi(str,perNams);  % prevent getting more than 1 value
                IColPer(ic)=ii(1);

                % Assert that this column exists
                if ~IColPer(ic)
                    error('Can''t find column <<%s>> in PER sheet',str);
                end
            end
        end

        % Populate cell array for all active stress periods.
        % First generate the list of cells to be used in every stress
        % period for this object.
        switch stressType
            case {'WEL','FLUX'}
                LIST = NaN(numel(o(io).Idx), 5);
            case {'GHB','DRN','CHD'}
                LIST = NaN(numel(o(io).Idx), 6);
            case 'RIV'
                LIST = NaN(numel(o(io).Idx), 7);
            otherwise
                % skip
        end
        LIST(:,2:4) = cellIndices([o(io).Idx],o(io).grSize,'LRC');

        % insert values from numerical columns, they are constant
        % between stress periods and therfore outside the loop
        % below
        for ic=1:numel(ICol)
            if ICol(ic)
                switch class(o)
                    case 'pointObj'
                        if Irel(ic)
                            LIST(:,4+ic) =  mf_interp2(gr.Xc,gr.Yc,gr.Z(:,:,1), ...
                                o(io).vertex(:,jx), o(io).vertex(:,jy)) ...
                                - o(io).vertex(:,ICol(ic)) ; % fill in the values
                        else
                            LIST(:,4+ic) =  o(io).vertex(:,ICol(ic)); % fill in the values
                        end
                    case 'lineObj'
                        if Irel(ic)
                            LIST(:,4+ic) = mf_interp2(gr.Xc,gr.Yc,gr.Z(:,:,1), ...
                                o(io).cellLineVals(:,jx), o(io).cellLineVals(:,jy)) ...
                                - o(io).cellLineVals(:,ICol(ic)); % fill in the values
                        else
                            LIST(:,4+ic) = o(io).cellLineVals(:,ICol(ic)); % fill in the values
                        end
                    case 'area2Obj'
                        if Irel(ic)
                            LIST(:,4+ic) = mf_interp2(gr.Xc,gr.Yc,gr.Z(:,:,1), ...
                                o(io).cellAreaVals(:,jx), o(io).cellAreaVals(:,jy)) ...
                                - o(io).cellAreaVals(:,ICol(ic)); % fill in the values
                        else
                            LIST(:,4+ic) = o(io).cellAreaVals(:,ICol(ic)); % fill in the values
                        end
                    otherwise
                        error('%s: illegal class <<%s>>',mfilename,class(o));
                end
            end
        end

        if strcmpi(o(io).name(1:5),'south')
            fprintf('\n');
        end
        % What to do with the extra array inputs ?
        % The adding to the relative drain elevation is obsolete because
        % this elevation is now dealt with in the spreadsheet.
        % TO 121011
        if ICol(1)
            if ismember(stressType,{'WEL','FLUX'})
                LIST(:,5) = LIST(:,5).*o(io).A(:).* o(io).V{1};
%                     else
%                         LIST(:,5) = LIST(:,5)             + o(io).V{1};
            end
        end
        if ICol(2)
            if strcmpi(stressType,'CHD')
%                        LIST(:,6) = LIST(:,6)             + o(io).V{2};
            else
                LIST(:,6) = LIST(:,6).*o(io).A(:).* o(io).V{2};
            end
        end
        if ICol(3)
%                   LIST(:,7) = LIST(:,7)                 + o(io).V{3};
        end

        % run over the stress periods to place the LIST and
        % get time varying values from the PER sheet

        hdrs = lower({'z','zBot','zRel','Hb','Db','H','D'});
        i = find( ismember( hdrs, lower(o(io).vertexHdr) ),1,'first');
        if isempty(i)
            error('Can''t find header matching any of { %s} in vertexHdr { %s}',...
                spintf(' %s',hdrs{:}),sprintf(' %s',o(io).vertexHdr));
        end
        iCol = strmatchi(hdrs{i},o(io).vertexHdr,'exact');
            
        switch class(o(io))
            case 'pointObj'
                switch i
                    case {1,2,4,6},  zObject =                   o(io).vertex(:,iCol);
                    case {5,7},      zObject = gr.Z(o(io).Idx) - o(io).vertex(:,iCol);
                    case 3,          zObject = gr.zRel2z(o(io).vertex(:,iCol),o(io).Idx);
                end
            case 'lineObj'
                zObject = [o(io).P.zm]';
                %zBottom = gr.ZBlay([o(io).P.idx]);
            case 'area2Obj'
                switch i
                    case {1,2,4,6}, zObject =                   o(io).cellAreaVals(:,iCol);
                    case {5,7},     zObject = gr.Z(o(io).Idx) - o(io).cellAreaVals(:,iCol);
                    case 3,         zObject = gr.zRel2z(o(io).cellAreaVals(:,iCol),o(io).Idx);
                end
        end

        if any(Irel)
%            idxTop = rem(o(io).Idx-1,gr.Nxy)+1;
            dem = gr.Z(gr.IdTop(o(io).Idx))';
        end
        
        for iper = numel(ACTIVE):-1:1
            
            if ~ACTIVE(iper), continue; end
            
            iCell = i+lastCell;
            LIST(:,1) = iper;

            % process values from char columns that refer to
            % columns in the PER sheet
            if IColPer(1),
                if ismember(stressType,{'WEL','FLUX'})
                    LIST(:,5) = perVals(iper,IColPer(1))*o(io).A(:).*o(io).V{1};
                else % elevation
                    if Irel(1)
                        LIST(:,5) = dem - perVals(iper,IColPer(1));
                    else
                        LIST(:,5) =       perVals(iper,IColPer(1));
                    end
                end
            end
            if IColPer(2)
                if strcmpi(stressType,'CHD')
                    if iRel(1)
                        LIST(:,6) = dem - perVals(iper,IColPer(2));
                    else
                        LIST(:,6) =       perVals(iper,IColPer(2));
                    end                        
                else
                    LIST(:,6) = perVals(iper,IColPer(2))*o(io).A(:).*o(io).V{2};
                end
            end
            if IColPer(3)  % Hb for RIV
                if Irel(3)
                   LIST(:,7) = dem - perVals(iper,IColPer(3));
                else
                   LIST(:,7) =       perVals(iper,IColPer(3));
                end
            end

            errMsg = [...
                '%s [partly?] above specified head in active stress period <<%d>>\n',...
                 'removing <<%d>> cells of a total of <<%d>>\n',...
                 'Verify graphcally what MODFLOW gets using function\n',...
                 '    gridObj/showWhatIsInMyStressFile(modflowStressFileName with extension .DRN)\n',...
                 'Make sure that elevations are > object or river bottom in the spreadsheet.\n',...
                 'Omitting zRel means drain is at ground surface!\n'];
            
            % Assert RIV stage is at least at riv bottom elevation
            % this is done in writeBCN, not here, to ensure it also works
            % with other RIV inputs.
            % TO 131119
            
            % TO 130927, to prevent specified drain elevation below
            % the actual elevation of te drain, use this
            if ismember(stressType,{'DRN' ,'RIV','CHD'})
                I = LIST(:,5)>=zObject;
                if all(~I)
                    warning(msgId,'%s ''%s'' of type %s in SP %d: head totally below object, entire object removed',...
                        class(o(io)),o(io).name,o(io).type,iper);
                    continue
                end
                if strcmp(stressType,'CHD')
                    I = I  & LIST(:,6)>=zObject;
                    if all(~I)
                        warning(msgId,'%s ''%s'' of type %s in SP %d: head totally below object, entire object removed',...
                            class(o(io)),o(io).name,o(io).type,iper);
                        continue
                    end
                end            
                if any(~I)
                    warning(msgId, errMsg,  o(io).type, o(io).name, iper, sum(~I), numel(I));
                    cellList{iCell,1}=LIST(I,:);
                else
                    cellList{iCell,1} = LIST;
                end
            else
                cellList{iCell,1} = LIST;
            end
            lastCell = numel(cellList);
        end
        if lastCell == 0
            if ~any(ACTIVE)
                error('No active stress periods exist for %s ''%s'' of type %s',...
                    class(o(io)),o(io).name,o(io).type);
            else
                error(['No cells to write for %s ''%s'' of type %s\n'...
                       'Remedy: see that all %s elevations are above the bottom\n',...
                       'of the cell in which they are placed and that the stress type\n',...
                       'has at least one active stress period'],...
                        class(o(io)),o(io).name,o(io).type,o(io).type);
            end
        end
        if strcmpi(stressType,'DRN'), warning('off',msgId); end
        warning('off',msgId);
    end
end
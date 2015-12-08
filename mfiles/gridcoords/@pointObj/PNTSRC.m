function cellList = PNTSRC(o,basename)
    %PNTSRC: generates PNTSRC for intput to MT3DMS and SEAWAT
    % cellList is one cell per stress period with the list of point
    % sources that have the stress period number in their first
    % column. The list as following structure on each line
    %  iPer L R C C1 ITYPE C1 C2 C3 ...
    %  where C is column of cell, and C1, C2, C3 are concentration
    %  of the species in the order of the CcolNames in the objects
    %  specified at the time of their creation. That is CcolNames
    %  are the names of the column headers in the input array of
    %  the object that hold the species. As always a numeric value
    %  implies a constant for all stress periods, whereas a string
    %  refers to the values per stress period in the column in the
    %  PER sheet that has a header equal to the value given in the
    %  input array in CcolNames.
    %
    % TO 130831

    [perNams,perVals,Nper] = getPeriods(basename);

    % Generates a cell array of length Nper, each of which has a list
    % of cells equal to length of all cells of all the objects
    % while their type translates into ITYPE
    cellList{Nper,1} = [];

    if isempty(o)
        return
    end            

    % Column numbers that correspond with CcolNames
    if ischar(o(1).CcolNames), o(1).CcolNames = {o(1).CcolNames}; end

    Nspec = numel(o(1).CcolNames);

    % Populate cell array for all, Nspec is number of chemical
    % species
    lastCell = 0;

    % object by object because objects may be active during
    % arbitrary stress periods
    for io=numel(o):-1:1

       % assert species column names are a cell array
       if ischar(o(io).CcolNames), o(io).CcolNames = {o(io).CcolNames}; end

        % see if a PER-sheet column was specified indicating when the
        % object is active.
       if ~isempty(o(io).activePerColName)
            if ismember(o(io).activePerColName(1),'~-')
                iActivePerCol = strmatchi(o(io).activePerColName(2:end),perNams);
            else
                iActivePerCol = strmatchi(o(io).activePerColName,perNams);
            end
            % Assert that this column exists in the PER sheet.
            if ~iActivePerCol
                error('can''t fine column <<%s>> in PER sheet to see when %s <<%s>> is active',...
                    o(io).activePerColName,o(io).class);
            end
        else
            iActivePerCol = 0;
       end

        for ic = Nspec:-1:1

            % species conc is given as a numeric value
           CCol_num(ic)=strmatchi(o(io).CcolNames{ic},o(io).vertexHdr,'exact');

           % species conc is given a string, referring to PER sheet column
           CCol_txt(ic)=strmatchi(o(io).CcolNames{ic},o(io).vertexTxtHdr,'exact');

           % Assert column existance in the PER sheet,
           % immediately get the corresponding column in the PER sheet
           if CCol_txt(ic)

               % get string from first line of o(io).vertexTxt
               perHdr = o(io).vertexTxt(1,CCol_txt(ic));

               % get the inex into of the corresponding column in the PER sheet
               CCol_txt(ic) = strmatchi(perHdr,perNams,'exact');

               % Assert it existance
               if ~CCol_text(ic)
                   error('Can''t fine column header <<%s>> in PER worksheet',perHdr);
               end
           end
        end

        % the list of SSM values for this object that will be used
        % in all stress periods that this object is active in.
        LIST = NaN(numel(o(io).Idx),4+2+Nspec);
        LIST(:,2:4) = cellIndices(o(io).Idx,o(io).grSize,'LRC');
        LIST(:,6  ) = eval(['itype.',o(io).type ';']);

        % Get the numeric values for the species involved
        for ic=1:Nspec
            if CCol_num(ic)  % column with numeric value
                switch class(o)
                    case 'pointObj'
                        LIST(:,6+ic)=o(io).vertex(:,CCol_num(ic));
                    case 'lineObj'
                        LIST(:,6+ic)=o(io).cellLineVals(:,CCol_num(ic));
                    case 'area2Obj'
                        LIST(:,6+ic)=o(io).cellAreaVals(:,CCol_num(ic));
                    otherwise
                        error('Illegal type <<%s>> object with name <<%s>> and type <<%s>>',...
                            class(o),o(io).name,o(io).type);
                end
            end
        end

       % get the stress periods in which this object is active
       if iActivePerCol
            if ismember(o(io).activePerColName(1),'~-')  % use opposite values
                IperActive = find(~logical(perVals(:,iActivePerCol)));
            else
                IperActive = find(logical(perVals(:,iActivePerCol)));
            end
       else
           IperActive = 1:Nper;
       end

       for i = numel(IperActive):-1:1

           iPer = IperActive(i);
           iCell= i+lastCell;

            cellList{iCell,1}      = LIST;
            cellList{iCell,1}(:,1) = iPer;

            % if species was defined as a string, then use the
            % corrresponding column from the PER sheet
            for ic=Nspec:-1:1
                iPerCol = CCol_txt(ic);
                if iPerCol
                    cellList{iCell}(:,6+ic)=perVals(i,iPerCol);
                end
            end

            % put species 1 also in column 5, this is for the case we
            % have only a single species (MT3DMS manual)
            cellList{iCell,1}(:,5)=cellList{iCell}(:,7);
            if Nspec==1
                cellList{iCell}=cellList{iCell}(:,1:6);
            end
       end
       lastCell = numel(cellList);
   end           
end

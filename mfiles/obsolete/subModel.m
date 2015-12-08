function [Model grNew] = subModel(dataFolder,xLim,yLim,zLim)
% restore this file from time machine

if iscell(dataFolder)
    Model = dataFolder; % dataFolder;
    
    for i=1:size(Model,1),  eval([ Model{i,1} '= Model{i,2};']); end

    GRID = Model{strmatchi('gridObj',Model(:,3)),2};

    if nargin<2
        xLim=GRID.xGr([1 end]);
        yLim=GRID.yGr([1 end]);
    end
    
    %% Indices to cutout the model
    Ix = find(GRID.xm>min(xLim) & GRID.xm<max(xLim)); %#ok<NASGU>
    Iy = find(GRID.ym>min(yLim) & GRID.ym<max(yLim)); %#ok<NASGU>

    if nargin<4,
        Iz=1:GRID.Nz; %#ok<NASGU>
    else
        Iz = find(GRID.zm>min(zLim) & GRID.zm<max(zLim)); %#ok<NASGU>
    end

    %% Here we use eval to load all the model arrays in a single cell array with fields
    % { ArrayVaraibleName   ActualArrayVariable ArrayVariableType }
    % The Array names are obtained from dir([dataloc '/*.mat'])
    % Each matfile has the name of the array it contains and has also the type
    % of the array on board. This tyype is also loaded into the Model cell
    % array. The type is used when handling the array when it comes to changing
    % an existing grid such as refining it.

    for i= 1:size(Model,1)
        %gr.cutout(VAR,VARtype,Ix,Iy)
        ie = num2str(i);
        fprintf(['Model{' ie ',2} = GRID.cutout(Model{' ie ',2},Model{',ie,',3},Ix,Iy,Iz);\n']);

        %% Variable itself and its type (3rd field)
        eval(   ['Model{' ie ',2} = GRID.cutout(Model{' ie ',2},Model{',ie,',3},Ix,Iy,Iz);']);
    end
    
else
    if isempty(strfind(dataFolder,'.mat'))
        error('%s: <<%s>> must end at ''.mat'' or ''/*.mat''',dataFolder);
    end

    %% Refer to a mat file containing all variables comprising the model + a gridObj
    if ~strfind( dataFolder,'*')
        %% Look for matfiles directory dataFolder

        load(dataFolder);

        varNames = who;

        if ~strmatchi('gridObj',varNames)
            error('The data mat file <<%s>> must contain a gridObj',dataFolder);
        end

        matdata(length(varNames),1).varName='Dummy';

        for i=1:length(varNames),
            matdata(i).varName = varNames{i};
            fprintf('Found variable <<%s>>\n',matdata(i).varName);
        end

    else  % wildcard given, so look for directory that contains all the variables in matfiles

        matdata = dir(dataFolder);

        if isempty(matdata)
            error('mfLab:subModel:fileNotFound',...
                '%s can''t find file(s) or <<%s>>',dataFolder);
        end

        bytes=0;

        dataFolder = regexp(dataFolder,'^[A-Za-z0-9./\\]+','match');

        for i=1:length(matdata)
            varName = regexp(matdata(i).name,'^[A-Za-z0-9]+','match');
            matdata(i).varName = varName{1};
            fprintf('Loading file %s ...',matdata(i).name);
            load(fullfile(dataFolder{1},matdata(i).name));    
            bytes = bytes + matdata(i).bytes;
            done
        end
        fprintf('Total number of bytes loaded = %ld = %.0f MB\n',bytes,bytes/1024/1024);
    end

    %% Evaluate the data, add them to the struct Model

    Model{length(matdata),3}='Dummy'; % Allocate
    
    %% Indices to cutout the model
    if ~exist('xLim','var')
        xLim = GRID.xGr([1 end]);
        yLim = GRID.yGr([1 end]);
    end
    Ix = find(GRID.xm>min(xLim) & GRID.xm<max(xLim)); %#ok<NASGU>
    Iy = find(GRID.ym>min(yLim) & GRID.ym<max(yLim)); %#ok<NASGU>
    if nargin<4,
        Iz=1:GRID.Nlay; %#ok<NASGU>
    else
        Iz = find(GRID.zm>min(zLim) & GRID.zm<max(zLim)); %#ok<NASGU>
    end

    %% Here we use eval to load all the model arrays in a single cell array with fields
    % { ArrayVaraibleName   ActualArrayVariable ArrayVariableType }
    % The Array names are obtained from dir([dataloc '/*.mat'])
    % Each matfile has the name of the array it contains and has also the type
    % of the array on board. This tyype is also loaded into the Model cell
    % array. The type is used when handling the array when it comes to changing
    % an existing grid such as refining it.

    for i= 1:length(matdata)
        % VAR = cutout(o,type,VAR,Ix,Iy,LAYCBDold,LAYCBDnew)
        ie = num2str(i);
        fprintf( ['[Model{' ie ',2}, Model{' ie ',3}]= GRID.cutout( ' matdata(i).varName ',' matdata(i).varName 'type ,Ix,Iy,Iz);' '\n']);

        %% Variable name
        Model{i,1}=matdata(i).varName;

        %% Variable itself and its type (3rd field)
        eval([ '[Model{' ie ',2}, Model{' ie ',3}]= GRID.cutout( ' matdata(i).varName ',' matdata(i).varName 'type ,Ix,Iy,Iz);']); 

    end
end

Model = updateWellSeries(Model);  % kind of cleanup

grNew = Model{strmatchi('gridObj',Model(:,3)),2};


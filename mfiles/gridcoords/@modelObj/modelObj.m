classdef modelObj
%MODELOBJ: generates objects to hold and manipulate model arrays and the grid that
%
% USAGE:
%      Model = modelObj(value,type)
%      Model = modelObj(ModelArray,xlim,ylim,Ilay)
%      Model = modelObj(dataFolder,xlim,ylim,Ilay)
%
% First form generates a modelObj for single array value of type 'type' 
% For instance
%     Model(10)=modelObj(IBOUND,'3Dlay');
%     Model(11)=modelObj(HK    ,'3Dlay');
%     Model(12)=modelOb( gr    ,'gridObj);
%     Model(13)=modelObj(wells ,'wellObj);
%     Model(14)=modelObj(DRN   ,'stress');
%     Model(15)=modelObj(PNTSRC,'stress');
%     Model(16)=modelObj(STCONC,3Dlay');
%     Model(17)=modelObj(RECH  ,3Dtime');
%
%  A model object array can hold the arrays and data arrays
%  that constitute a model as defined in mf_adapt. Holding the
%  model in an array of modelObjects allows splitting and merging
%  layers, including stress arrays (like WEL, DRN, RIV, GHB, CHD)
%  and the grid.
%  Each modelObj has 2 fields: value and type.
%  The value is a variablelike HK, gr, DRN and its type indicates
%  what type of model array it is.
%
%  Type can be oneof
%   'zlist'      a variable with one value per layer (see LAY worksheet)','WETDRY, LAYCON, ...';
%   '3Dlay'      a 3D array of size (Ny,Nx,Nlay)','HK, STRTHD, ...';
%   '3Dcbd'      a 3D array of size (Ny,Nx,Ncbd)','VCONT, VKCB ...';
%   '3Dtime'     a 3D time array of size (Ny,Nx,NPER)','RECH, EVAP, ...';
%   'stress'     a list array for a stress','DRN, GHB, RIV, CHD, STR, ...';
%   'gridObj'    a gridObj
%   'wellObj'       array of wellObj
%   'wellSeriesObj' array with wellSeriesObj
%   'kDObj'         array with kD objects (constructed from HNI well file)
%   'modelObj'   a complete array of modelObj
%   'struct'     a struct like from readDat,readBud,readMt3d
%
%
%  The tird form generates a submodel from a larger model that is contained
%  dataFolder (stored in mat files) in the data directory. Each of the mat
%  file names in that directory corresponds to the variables the mat file
%  contains (one variable per matfile).
%  The directory must include a gridObj named gr or GRID.
%
%  xlim is the part of the x-axis to select (xmin xmax) of new grid.
%  ylim is the part of the y-axis to select (ymin ymax) of new grid
%  Ilay   are the layers of the new grid (subection of those of
%  the old grid.
%
% Contents of the mat files:
% Each mat file contains two parameters,
%  1) the Matlab variable
%  2) its mfLab type (oneof the list shown above)
%
%
%% Converting Block Centered Flow Package (BCF) dataFolder to Layer Property Flow Package (LFP) data.
% If we intend transport modeling, we must use the LPF package instead of BCF.
% Moreover, the use of the confining bed option in the LPF
% package is discouraged with MT3DMS (Zheng, 2005).
% Therefore, we will convert confining beds into regular model layers.
% Therefore, when we cutout a submodel, BCF arrays (TRAN, HY, VCONT)
% will be transferred to LPF arrays (HK, VK, VKCB).
%
% Note that mfLab prefers the use of VK for vertical cell conductance instead
% of VKA as spcecified in the MODFLOW 2000 LPF manual. This is for user
% convenience and to avoid the confusion that may go along with the use
% of VKA (combined with LAYVKA) to decide on a layer by layer basis whether
% VKA is the vertical conductivity or the vertical anisotropy. However,
% the use of VKA combined with VK is posissible with mfLab.
% Whenever mfLab encounters VKA in the workspace it checks LAYVKA in the
% workbook LAY sheet to decide between vertical conductivity and vertical
% anisotropy. If, however, mfLab encounters VK, it does not check LAYVKA
% but just always uses VK as vertical hydraulic conductivity.
%
% Conversion one on grid object read from the data directory.
% The conversion requires the grid. This grid must be in one of the mat
% files in the data directory. The name of this mat file is either gr.mat
% or GRID.mat and must have the associated type 'gridObj' stored along with the
% gridObj in the same mat file.
%
% A resulting submodel will be stored in the workspace variable of class modelObjModel
% in the file  Model.mat
%
% TO 120607 120810
            

    properties
        dataFolder,        % origin of the data
        description, % description of model
        name,        % varname
        var,         % variable (array)
        type,        % mfLab type
        UserData     % to be used by user for his own purposes
    end
    methods
        function o = modelObj(varargin)
            %MODELOBJ: constructor of modelObj
            %
            % USAGE:
            %      Model = modelObj(value,type)
            %      Model = modelObj(ModelArray,xlim,ylim,Ilay)
            %      Model = modelObj(dataFolder,xlim,ylim,Ilay)
            %
            % TO 120810
            
            if nargin == 0
                % to generate N objects use
                % Model(N) = modelObj();
                % which is a model array which has N variables.
                return;
            end
            
            %% call:modelObj(var,type)
            if ~isempty(inputname(1)) && length(varargin)>1 && ischar(varargin{2})
                for io=numel(o):-1:1
                    o(io).name = inputname(1);
                    o(io).var  = varargin{1}; % ylim now used for var
                    o(io).type = varargin{2}; % ylim now used for type

                    if ~strmatchi(o(io).type,{'zlist','3Dlay','3Dcbd','3Dtime','stress',...
                            'gridObj','wellObj','wellSeriesObj','kDObj','struct'});
                        error('illegal type <<%s>>>',o.type);
                    end
                end
                return;
            end
            
            
            %% if input is a struct with fields as in modelObj
            % then turn it into an array of modelObj
            if isstruct(varargin{1})
                b = varargin{1};
                if all(ismember({'name','var','type'},fieldnames(b)))
                    for io=numel(b):-1:1
                        o(io)=modelObj(b(io).name,b(io).var,b(io).type);
                    end
                end
                return;
            end
                        
            % At this point we expect a directory pathname as first argument
            % possibly followed by xLim,yLim,ILay.
            % The data should be in that folder in matfiles.
            % The matfiles that will be read and interpreted as individual
            % moddel arrays.
            % Each mat file should contain a variable and its type
            % (i.e. '3Dlay', '3Dtime' etc according to the list given above.
            [dataFolder_,varargin] = getNext(varargin,'char','');
            if ~exist(dataFolder_,'dir')
                error('At this point, varargin{1}=<< %s >> must be a path/file name',dataFolder_);
            end
                
            if numel(varargin)>=1, xlim = varargin{1}; else xlim=[-Inf,Inf]; end
            if numel(varargin)>=2, ylim = varargin{2}; else ylim=[-Inf,Inf]; end

            %% Check for the gridObj in file gr*.mat (case insensitive)

            d = dir(fullfile(dataFolder_,'*.mat'));
            
            % Skip the WEL.mat file (WEL lists) because we have well.mat
            % (well objects)
            % Skip NHIvars because it does not contain an array
            d(ismember({d.name},{'WEL.mat','NHIvars.mat'}))=[];

            % put grid in front
            iGr = strmatchi('gr.mat'  ,{d.name},'exact');
            if iGr
                dg = d(iGr);
                d(iGr)=[];
                d = [dg; d];
                iGr = strmatchi('GRID.mat',{d.name},'exact');
                if iGr
                    d(iGr)=[];
                end
            else
                iGr = strmatchi('GRID.mat',{d.name},'exact');
                if iGr
                    dg = d(iGr);
                    d(iGr)=[];
                    d = [dg; d];
                else
                    error('Found no matfile which contains the grid in <<%s>>',dataFolder_);
                end
            end

            % load the grid, is now first in row
            load([dataFolder_ d(1).name]);
            if exist('GRID','var'), gr=GRID; clear('GRID'); end

            % having the grid compute the indices of the selected model
            Ix = between(gr.xm,xlim);
            Iy = between(gr.ym,ylim);

            if nargin<4
                Ilay=1:gr.Nlay;
            else
                Ilay = varargin{4};
                Ilay = Ilay(ismember(Ilay,1:gr.Nz));
            end
                
            % grid will be cutOut last
            for i=length(d):-1:1
                load([dataFolder_ d(i).name]);
                [~,varName] = fileparts(d(i).name);

                % Verify existance of var and varType
                varType = [varName 'type'];
                if ~exist(varName,'var')
                    error('%s: file <<%s>> does not contain variable <<%s>>',mfilename,d(i).name,varName);
                end
                if ~exist(varType,'var')
                    error('%s: file <<%s>> does not contain variable <<%s>>',mfilename,d(i).name,varType);
                end

                % Show what was found
                fprintf('Model.varName = %s,  Model.varType = %s\n',varName,varType);

                % Cutout and store the submodel for this variable
                o(i).dataFolder = dataFolder_;
                o(i).name = varName;
                o(i).var  = eval(varName); 
                o(i).type = eval(varType);
                o(i)      = gr.cutout(o(i),Ix,Iy,Ilay);

                if ~strcmp(o(i).type,'gridObj')
                    clear(varName,varType);
                end
                
            end
            o(1).name = 'gr';
        end
        function gr = grid(o)
            % gr = Model.grid() -- issue the grid
            i = strmatchi('gridObj',{o.type});
            if ~i || numel(i)>1
                error('gridObj in <<%s>> ambiguous',o.type);
            else
                gr = o(i).var;
            end
        end
        function o = descr(o,txt)
            % Model = Model.descr(txt) - store txt in model as a
            % description
            if nargin<2
                fprintf('Actions taken to generate this model:\n');
                for i=1:size(o(1).description,1);
                    fprintf('%2d: %s\n',i,o(1).description{i});
                end
            else
                for i=1:length(o)
                    if isempty(o(i).description)
                        o(i).description = {txt};
                    else
                        o(i).description{end+1,1} = txt;
                    end
                end
            end
        end        
    end
    methods (Static)
        function tp = types()
            % Model.types() -- prints the legal mfLab variable types in the
            % Model object.
            
             tp = {
            'zlist'        ,'variable with one value per layer (see LAY worksheet)','WETDRY, LAYCON, ...';
            '3Dlay'        ,'3D array of size (Ny,Nx,Nlay)','HK, STRTHD, ...';
            '3Dcbd'        ,'3D array of size (Ny,Nx,Ncbd)','VCONT, VKCB ...';
            '3Dtime'       ,'3D time array of size (Ny,Nx,NPER)','RECH, EVAP, ...';
            'stress'       ,'list array for a stress','DRN, GHB, RIV, CHD, STR, ...';
            'gridObj'      ,'grid object','gr or GRID';
            'wellObj'      ,'array with well objects','well ...';
            'wellSeriesObj','array with wellSeriesObj','wellSeries ...';
            'kDObj'        ,'array with kD objects (constructed from HNI well file)','kD';
            'modelObj'     ,'a complete model objects','Model1 Model2 etc';
            'struct'       ,'a struct like from readDat,readBud,readMt3d','*.hds, *.ddn, *.bgt, MT3D00?.UCN,...'};

            fprintf('\nLegal mfLab types associated with the variables in the mat files:\n');
            fprintf('%-13s     %-55s  %s\n','mfLabType','description','example variables of this type');
            for i=1:size(tp,1)
                fprintf('%-13s     %-55s, %s\n',tp{i,1:3});
            end
            fprintf('\n');
        end
    end
end
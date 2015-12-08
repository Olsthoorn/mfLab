classdef resObj
	%res = resObj(basename,grid)
    %	generates 4 commonly used functions of mf_analyze.
    
    %   Streamlines, plots the streamlines/flowlines, based on the
    %   flowrightface term. Lines of equal concentration, plots planes with 
    %   equal concentration values. Concentrationslices, uses buildin slice
    %   function to plot a slice on a certain axis. contourf, is the 
    %   buildin contourf-function, but slighlty changed so it can be
    %   plotted in a 3d enviroment.
    %   All plots are being send to the current axis
    
    %   If needed, the C, H, B values can be used directly.
    
    %   As for plot arguments, all patch-plot-arguments are supported.   
      
    %	todo:
    %   - combine the seperate direction slices to more general functions.
    %   one function for either X,Y and Z
    %	- flow direction arrows
    %	- more efficient updating of plots, instead of deleting handle.
    %	Although this is already fast enough.
    %	- working with multiple concentrations
    %	- specify axes, as a variable in the entire object
    %   - make movie function
    %   - certain properties need to be set to hidden. eg. plot handles
    %	120623, 120724, 120801, 120911  b.f.destombe@student.tudelft.nl
    
    properties (Constant)
        type='resObj';
    end
    properties (Dependent=false,SetAccess=public)
        cRange=[0 1];           % needs to be set to setup
        time=1;                 % current timestep for all plots, first timestep is defined, use update function for other timesteps
                
        strmRange               % these are the values for each flowline
        strmPlotarg={};         % containing the plot vargins
                
        equiConcRange           % these are the values for each flowline
        equiConcPlotarg={};     % containing the plot vargins
        
        concSliceX=[];          % a concentration slice/plane on x-ax
        concSliceY=[];          % a concentration slice/plane on y-ax
        concSliceZ=[];          % a concentration slice/plane on z-ax
        concSlicePlotarg={};    % containing the plot vargins
        
        concfX=[];              % a concentration slice/plane on x-ax, using contourf
        concfY=[];              % a concentration slice/plane on y-ax, using contourf
        concfZ=[];              % a concentration slice/plane on z-ax, using contourf
        concfArg=8;             % set contourf option; n or v in the contourf helpfile
        concfPlotarg={};        % containing the plot vargins
    end
    properties (Dependent=false,SetAccess=private)
        basename, gr
        
        bMask =0.001;           % to prevent roundoff errors in the budgetfile
        bThreshold= 0.999;
        
        strmPlot                % containing the handles
        equiConcPlot            % containing the handles
        concSlicePlot           % containing the handles
        concfPlot               % containing the handles
    end
    properties (Dependent=true, SetAccess=private)
        C, H, B                 % Read datafiles
        Bright, Bfront, Blower  % Read budgetfile with different flowfaces
        U                       % custom meshgrid of the cellnodes
        
        perData, nper, perlen   % get period data, excel is only read once. If needed it is easy to add your own variables, use perlen as template
    end
    methods
        function o = resObj(basename,grid)
            if nargin ~= 2
                error('not the right input arguments for resObj: use resObj(basename,grid)');
            end
            
            if isempty(basename),   display(basename);  error('empty basename, check it in mf_adapt.\n');end
            if isempty(grid),       display(grid);      error('empty gridObj, check it in mf_adapt.\n'); end
            
            if ~ischar(basename),   display(basename);  error('basename is not a string');  end
            if ~isa(grid,'gridObj'),display(grid);      error('the grid is not a gridObj'); end
                        
            o.basename  = basename;
            o.gr        = grid;
        end

        function C = get.C(o)
            persistent saved_C
            if isempty(saved_C)
                saved_C   = readMT3D('MT3D001.UCN');
                saved_C   = maskHC(saved_C,[0 1],[0 1]);
                for it  = 1:length(saved_C)   % convert concentrationvalues to added-salt-content [kg/m3]
                    saved_C(it).values     = saved_C(it).values.*(o.cRange(end)-o.cRange(1))+o.cRange(1);
                end
            end
            C   = saved_C;
        end
        
        function H = get.H(o)
            persistent saved_H
            if isempty(saved_H)
                saved_H   = readDat([o.basename,'.hds']);
            end
            H   = saved_H;
        end
        
        function B = get.B(o)
            persistent saved_B
            if isempty(saved_B)
                saved_B  = mf_Psi(readBud([o.basename,'.BGT']));
                for it=1:length(saved_B)
                    for iLbl=1:length(saved_B(it).label)
                        saved_B(it).term{iLbl}(abs(saved_B(it).term{iLbl})<o.bThreshold*o.bMask)=0;
                    end
                end
            end
            B   = saved_B;
        end
        
        function Bright = get.Bright(o)
            persistent saved_Bright
            if isempty(saved_Bright)
                saved_Bright  = mf_Psi(readBud([o.basename,'.BGT'],{'FLOWRIGHTFACE'}));
                for it=1:length(saved_Bright)
                    for iLbl=1:length(saved_Bright(it).label)
                        saved_Bright(it).term{iLbl}(abs(saved_Bright(it).term{iLbl})<o.bThreshold*o.bMask)=0;
                    end
                end
            end
            Bright   = saved_Bright;
        end
        function Bfront = get.Bfront(o)
            persistent saved_Bfront
            if isempty(saved_Bfront)
                saved_Bfront  = mf_Psi(readBud([o.basename,'.BGT'],{'FLOWFRONTFACE'}));
                for it=1:length(saved_Bfront)
                    for iLbl=1:length(saved_Bfront(it).label)
                        saved_Bfront(it).term{iLbl}(abs(saved_Bfront(it).term{iLbl})<o.bThreshold*o.bMask)=0;
                    end
                end
            end
            Bfront   = saved_Bfront;
        end
        function Blower = get.Blower(o)
            persistent saved_Blower
            if isempty(saved_Blower)
                saved_Blower  = mf_Psi(readBud([o.basename,'.BGT'],{'FLOWLOWERFACE'}));
                for it=1:length(saved_Blower)
                    for iLbl=1:length(saved_Blower(it).label)
                        saved_Blower(it).term{iLbl}(abs(saved_Blower(it).term{iLbl})<o.bThreshold*o.bMask)=0;
                    end
                end
            end
            Blower   = saved_Blower;
        end
         
        function U = get.U(o), [U.X,U.Y,U.Z]=meshgrid(o.gr.xc,o.gr.yc,o.gr.zc); end
                
        function perData = get.perData(o)
            persistent saved_perData
            if isempty(saved_perData)
                [pernams,pervals] = getPeriods(o.basename);
                saved_perData = cell2struct(num2cell(pervals'),genvarname(pernams),1);
            end
            perData = saved_perData;
        end
        
        function nper = get.nper(o)
            nper = o.perData(1).nPer;
        end
        
        function perlen = get.perlen(o)
        % Use this function as template if you'd need other data from the
        % PER-excelsheet
            perlen=zeros(1,length(o.perData));
            for i=1:length(o.perData)
                perlen(i)=o.perData(i).PERLEN;
            end
        end
        
        function o = set.cRange(o,range)
            if nargin ~= 2,
                error([ 'not the right input arguments for setCrange: use setCrange(cRange).\n',...
                    'example1: cRange=[0 29]\n',...
                    'example2: cRange=[29]; This assumes crange is from 0 to 29\n']);
            end
            if ~isa(range,'double'), error('range is not a double'); end
            switch numel(range)
                case 1
                    o.cRange = unique([0 range]);
                case 2
                    o.cRange = unique(range);
                otherwise
                    error(['wrong size of cRange\n',...
                        'example1: cRange=[0 29]\n',...
                        'example2: cRange=[29]; This assumes crange is from 0 to 29\n']);
            end
%             set(gca,'cLim',o.cRange);
        end
        
        %% Plot functions
        %   general plot properties
        function o = set.time(o,timestep)
            %set the period you want to plot.
            %   when plot method, ex. strmPlot, are called directly, this
            %   timestep will be plotted.
            %   timestep needs to be between 1 and resObj.nper. remember
            %   that the update method can include a timestep, which makes
            %   it unnessesarry to set this manually.
            %   Example: resObj.time(5);
            if ~isa(timestep,'double'), error('timestep is not a double'); end
            if numel(timestep)~=1, error('please enter just one timestep');end
            o.time = unique(timestep);
        end
        
        function lim(o)
            %sets x,y,z lim from gridobj to gca;loads c-lim from cRange
            %   Use this method in the beginning of your mf_analyze,
            %   directly after defining resObj.cRange.
            %   Example: resObj.lim;
            xlabel('x')
            ylabel('y')
            zlabel('z')
            
            set(gca,'YLim',o.gr.ylim)
            set(gca,'XLim',o.gr.xlim)
            set(gca,'ZLim',o.gr.zlim)
            set(gca,'CLim',o.cRange)
            
            set(gca,'NextPlot','add')
            view([0 -1 0])
        end
        function update(o,varargin)
            %plot all set plotcommands.
            %   usage1: resObj.update;
            %           plots all for which a range/slice location is set
            %           for the period set in resObj.time
            %   usage2: resObj.update(5);
            %           plots all for which a range/slice location is set
            %           for period 5.
            
            switch nargin-1
                case 0
                case 1
                    timestep=varargin{1};
                    if ~isa(timestep,'double'), error('timestep is not a double'); end
                    if numel(timestep)~=1, error('please enter just one timestep');end
                    o.time = unique(timestep);
                otherwise
                    error('please use update function with no arguments or one argument, the timestep');
            end            
            o.strmUpdate;
            o.equiConcUpdate;
            o.concSliceUpdate;
            o.concfUpdate;
        end
        
        %  streamlines
        function o = set.strmRange(o,range)
            if ~isa(range,'double'), error('range is not a double'); end
%             if isempty(range), error('equiConcRange is empty, nothing to plot'); end
            o.strmRange = unique(range);
        end
        function o = set.strmPlotarg(o,varargin)
            %these arguments are being passed to the patchobjects.
            %   all patch-plot-arguments are supported.
            %   Example: resObj.strmPlotarg={'FaceAlpha',0.5,'EdgeAlpha',0.2};
            
            if ~iscell(varargin), error('arguments are not in cellformat'); end
            o.strmPlotarg = varargin{:};
        end
        
        function strmUpdate(o,varargin)
            switch nargin-1
                case 0
                case 1
                    timestep=varargin{1};
                    if ~isa(timestep,'double'), error('timestep is not a double'); end
                    if numel(timestep)~=1, error('please enter just one timestep');end
                    o.time = unique(timestep);
                otherwise
                    error('please use update function with no arguments or one argument, the timestep');
            end
            delete(findobj('DisplayName','streamlines'));
            if ~isempty(o.strmRange),     o.strmPlot;     end
        end
        function strmPlot = get.strmPlot(o)
            if ~isempty(o.strmRange)
                hold on
                V   = cumsum(o.Bright(o.time).term{1},3); % take the cumsum over the z-ax
                strmPlot = NaN(1,length(o.strmRange));
                for i = 1:length(o.strmRange)
                    strmPlot(i) = patch(isosurface(o.U.X,o.U.Y,o.U.Z,V,o.strmRange(i)));
                    isonormals(o.U.X,o.U.Y,o.U.Z,V,strmPlot(i));
                    set(strmPlot(i),'FaceColor','interp','EdgeColor','interp',o.strmPlotarg{:});
                    set(strmPlot(i),'DisplayName','streamlines');
                    isocolors(o.U.X,o.U.Y,o.U.Z,o.C(o.time).values,strmPlot(i));
                end
                drawnow;
            else
                warning('resObj:resObj:emptystrmRange','please define strmRange');
                strmPlot = [];
            end
        end
              
        %   planes of equal concentration
        function o = set.equiConcRange(o,range)
            %for each given concentration, a plane will be drawn for each
            %given concentration.
            %   The equiconcrange needs to be within cRange
            %   Example: resObj.equiConcRange=0:.1:1;
            
            if ~isa(range,'double'), error('range is not a double'); end
%             if isempty(range), error('equiConcRange is empty, nothing to plot'); end
            o.equiConcRange = unique(range);
        end        
        function o = set.equiConcPlotarg(o,varargin)
            %these arguments are being passed to the patchobjects.
            %   all patch-plot-arguments are supported.
            %   Example: resObj.equiConcPlotarg={'FaceAlpha',0.5,'EdgeAlpha',0.2};
            
            if ~iscell(varargin), error('arguments are not in cellformat'); end
            o.equiConcPlotarg = varargin{:};
        end
        
        function equiConcUpdate(o,varargin)
            %plot all planes of equal concentraion set in
            %resObj.equiConcRange.
            %   usage1: resObj.equiConcUpdate;
            %           plots the period set in resObj.
            switch nargin-1
                case 0
                case 1
                    timestep=varargin{1};
                    if ~isa(timestep,'double'), error('timestep is not a double'); end
                    if numel(timestep)~=1, error('please enter just one timestep');end
                    o.time = unique(timestep);
                otherwise
                    error('please use update function with no arguments or one argument, the timestep');
            end
            delete(findobj('DisplayName','lines of equal concentration'));
            if  ~isempty(o.equiConcRange),  o.equiConcPlot;             end
        end        
        function equiConcPlot = get.equiConcPlot(o)
            if ~isempty(o.equiConcRange)
                if any((o.equiConcRange < o.cRange(1) | o.equiConcRange > o.cRange(2)))
                    error('The equiConcRange values need to be within resObj.cRange.\n so:higher then %f and lower then %f. Or set resObj.cRange',o.cRange(1),o.cRange(2))
                end
                hold on
                equiConcPlot = NaN(1,length(o.equiConcRange));
                for i = 1:length(o.equiConcRange)
                    equiConcPlot(i) = patch(isosurface(o.U.X,o.U.Y,o.U.Z,o.C(o.time).values,o.equiConcRange(i)));
                    isonormals(o.U.X,o.U.Y,o.U.Z,o.C(o.time).values,equiConcPlot(i));
                    set(equiConcPlot(i),'FaceColor','interp','EdgeColor','interp',o.equiConcPlotarg{:});
                    set(equiConcPlot(i),'DisplayName','lines of equal concentration');
                    isocolors(o.U.X,o.U.Y,o.U.Z,o.C(o.time).values,equiConcPlot(i))
                end
                drawnow;
            else
                warning('resObj:resObj:emptyequiConcRange','please define equiConcRange');
                equiConcPlot = [];
            end
        end
        
        %   concentrationslices
        function o = set.concSliceX(o,sliceX)
            if ~isa(sliceX,'double'), error('X is not a double'); end
            o.concSliceX = unique(sliceX);
        end
        function o = set.concSliceY(o,sliceY)
            if ~isa(sliceY,'double'), error('Y is not a double'); end
            o.concSliceY = unique(sliceY);
        end
        function o = set.concSliceZ(o,sliceZ)
            if ~isa(sliceZ,'double'), error('Z is not a double'); end
            o.concSliceZ = unique(sliceZ);
        end
        function o = set.concSlicePlotarg(o,varargin)
            %these arguments are being passed to the patchobjects.
            %   all patch-plot-arguments are supported.
            %   Example: resObj.concSlicePlotarg={'FaceAlpha',0.5,'EdgeAlpha',0.2};
            
            if ~iscell(varargin), error('arguments are not in cellformat'); end
            o.concSlicePlotarg = varargin{:};
        end
        
        function concSliceUpdate(o,varargin)
            switch nargin-1
                case 0
                case 1
                    timestep=varargin{1};
                    if ~isa(timestep,'double'), error('timestep is not a double'); end
                    if numel(timestep)~=1, error('please enter just one timestep');end
                    o.time = unique(timestep);
                otherwise
                    error('please use update function with no arguments or one argument, the timestep');
            end
            delete(findobj('DisplayName','concentrationslices'));
            if (~isempty(o.concSliceX)||~isempty(o.concSliceY)||~isempty(o.concSliceZ)), o.concSlicePlot; end
        end        
        function concSlicePlot = get.concSlicePlot(o)
            if (~isempty(o.concSliceX)||~isempty(o.concSliceY)||~isempty(o.concSliceZ))
                if any(o.concSliceX < min(o.gr.xc))
                    error('the resObj.concSliceX slice value needs to be higher then %f and lower then %f.',min(o.gr.xc),max(o.gr.xc))
                end
                if any(o.concSliceY < min(o.gr.yc))
                    error('the resObj.concSliceY slice value needs to be higher then %f and lower then %f.',min(o.gr.yc),max(o.gr.yc))
                end
                if any(o.concSliceZ < min(o.gr.zc))
                    error('the resObj.concSliceZ slice value needs to be higher then %f and lower then %f.',min(o.gr.zc),max(o.gr.zc))
                end
                hold on
                concSlicePlot = slice(o.U.X,o.U.Y,o.U.Z,o.C(o.time).values,o.concSliceX,o.concSliceY,o.concSliceZ);
                for i=1:length(concSlicePlot)
                    set(concSlicePlot(i),'FaceColor','interp','EdgeColor','interp','DiffuseStrength',1,o.concSlicePlotarg{:});
                    set(concSlicePlot(i),'DisplayName','concentrationslices');
                end
                drawnow;
            else
                warning('resObj:resObj:emptyconcSlice','please define concSliceX, concSliceY, concSliceZ to your likings');
                concSlicePlot = [];
            end
        end
        %   contourf3dslices
        function o = set.concfX(o,sliceX)
            if ~isa(sliceX,'double'), error('X is not a double'); end
            o.concfX = unique(sliceX);
        end
        function o = set.concfY(o,sliceY)
            if ~isa(sliceY,'double'), error('Y is not a double'); end
            o.concfY = unique(sliceY);
        end
        function o = set.concfZ(o,sliceZ)
            if ~isa(sliceZ,'double'), error('Z is not a double'); end
            o.concfZ = unique(sliceZ);
        end
        function o = set.concfArg(o,v)
            %	contourf(Z,n) draws a filled contour plot of matrix Z with 
            %   n contour levels.
            %	contourf(Z,v) draws a filled contour plot of matrix Z with 
            %   contour lines at the data values specified in the 
            %   monotonically increasing vector v. The number of contour 
            %   levels is equal to length(v). To draw a single contour of 
            %   level i, use contour(Z,[i i])
            if ~isa(v,'double'), error('concfArg is not a double'); end
            o.concfArg = v;
        end            
        function o = set.concfPlotarg(o,varargin)
            %these arguments are being passed to the patchobjects.
            %   all patch-plot-arguments are supported.
            %   Example: resObj.concfPlotarg={'FaceAlpha',0.5,'EdgeAlpha',0.2};
            
            if ~iscell(varargin), error('arguments are not in cellformat'); end
            o.concfPlotarg = varargin{:};
        end
        
        function concfUpdate(o,varargin)
            switch nargin-1
                case 0
                case 1
                    timestep=varargin{1};
                    if ~isa(timestep,'double'), error('timestep is not a double'); end
                    if numel(timestep)~=1, error('please enter just one timestep');end
                    o.time = unique(timestep);
                otherwise
                    error('please use update function with no arguments or one argument, the timestep');
            end
            delete(findobj('DisplayName','contourf3dslices'));
            if (~isempty(o.concfX)||~isempty(o.concfY)||~isempty(o.concfZ)), o.concfPlot; end
        end
        function concfPlot = get.concfPlot(o)
            % becareful, plots the concentrations of the nearest cellnodes.
            % it does not interpolates inbetween cellnodes. future option. 
            if (~isempty(o.concfX)||~isempty(o.concfY)||~isempty(o.concfZ))
                if isempty(o.concfArg),error('resObj.concfArg is not supposed to be empty');end
                hold on
                concfPlot = [];
                if ~isempty(o.concfX)
                    cax=gca; % define current axis, for future use
                    for i=1:length(o.concfX)
                        if (o.concfX<min(o.gr.xc) || o.concfX>max(o.gr.xc))
                            error('the resObj.concfX slice value needs to be higher then %f and lower then %f.',min(o.gr.xc),max(o.gr.xc))
                        end
                        x=o.concfX(i);
                        y=o.gr.yc;
                        z=squeeze(o.gr.zc);
                        xi=round(interp1(o.gr.xc,1:length(o.gr.xc),x));
                        c=squeeze(o.C(o.time).values(:,xi,:));
                        hand=contourf3dX(cax,x,y,z,c,o.concfArg,o.concfPlotarg{:});
                        for ii = 1:length(hand)
                            set(hand(ii),'DisplayName','contourf3dslices');
                        end
                        concfPlot = [concfPlot;hand]; %#ok<AGROW>
                    end
                end
                if ~isempty(o.concfY)
                    cax=gca; % define current axis, for future use
                    for i=1:length(o.concfY)
                        if (o.concfY<min(o.gr.yc) || o.concfY>max(o.gr.yc))
                            error('the resObj.concfY slice value needs to be higher then %f and lower then %f.',min(o.gr.yc),max(o.gr.yc))
                        end
                        x=o.gr.xc;
                        y=o.concfY(i);
                        z=squeeze(o.gr.zc);
                        yi=round(interp1(o.gr.yc,1:length(o.gr.yc),y));
                        c=squeeze(o.C(o.time).values(yi,:,:));
                        hand=contourf3dY(cax,x,y,z,c,o.concfArg,o.concfPlotarg{:});
                        for ii = 1:length(hand)
                            set(hand(ii),'DisplayName','contourf3dslices');
                        end
                        concfPlot = [concfPlot;hand]; %#ok<AGROW>
                    end
                end
                if ~isempty(o.concfZ)
                    cax=gca; % define current axis, for future use
                    for i=1:length(o.concfZ)
                        if (o.concfZ<min(o.gr.zc) || o.concfZ>max(o.gr.zc))
                            error('the resObj.concfZ slice value needs to be higher then %f and lower then %f.',min(o.gr.zc),max(o.gr.zc))
                        end
                        x=o.gr.xc;
                        y=o.gr.yc;
                        z=o.concfZ(i);
                        zi=round(interp1(o.gr.zc,1:length(o.gr.zc),z));                        
                        c=squeeze(o.C(o.time).values(:,:,zi));
                        hand=contourf3dZ(cax,x,y,z,c,o.concfArg,o.concfPlotarg{:});
                        for ii = 1:length(hand)
                            set(hand(ii),'DisplayName','contourf3dslices');
                        end
                        concfPlot = [concfPlot;hand]; %#ok<AGROW>
                    end
                end
                drawnow;
            else
                warning('resObj:resObj:emptyconcf','please define concfX, concfY, concfZ to your likings');
                concfPlot = [];
            end
        end
    end
end
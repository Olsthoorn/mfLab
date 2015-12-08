classdef animateObj
%ANIMATEOBJ class def for animate obj
%
%  USAGE:
%     first request an animateObj, then visualize using one of its methods.
%     animate = animateObj(basename,gr,titleStrs,'mask',patches);
%  INPUT
%     basename as usual
%     gr gridObj
%     titleStrs cell array of titles to tell which data to load for subsequent
%          recognized 'head',drawd','budget',.
%          other are names of the concentrations in MT3D00?.UCN like
%            'Salinity','Tracer'
%     'mask',patches set of patches to mask the drawing, where patches are
%     a struct array with fields 'x','y','color'.
%     the order follows the order in which species nrs are used in Seawat.
% OUTPUT:
%     animateObj
% EXAMPLE:
%     animate =
%     animateObj(basename,gr,{'temp','salinity','heads',budget'},'mask,patches);
%            'temp' and 'salinity' are the names of the 'concentrations' in
%            MT3D001.UCN and MT3D002.UCN respectively.
% Notices:
%     type the name of the animateObj to see its methods
%
% TO 120420 120901
    properties
        gr                   = []; % the grid
        NCOMP                = 1;
        STCONC               = {[]};
        basename             ='Unknown!';
        titleStr             ='Title string not initialized!';
        numberOfHeadContours = 50;
        numberOfConcContours = 50;
        numberOfStreamLines  = 50;
        xLim
        yLim
        framerate            = 15;
        quality              = 80;
        time      % time. We can always add 
        H         % heads
        D         % drawdown
        UCN       % concentrations or temperatures from readDat or readMT2D
        % time of simulation, end of each stress period. Add t0 as datenum
        % to obtain datenumes here
        t0   = 0; % start time as datenum or zero if not used
        B         % budget
        %well     % wells
        tscale = 1;   % time conversion factor in titles
        tdim   = 'd'; % time dimension string in titles
        patch  = [];  % a mask/patch containing patch.x, patch.y, patch.color
        psiMask= [];  % to mask the stream function (fields x,y,)
        dateformat = 'dd/mmm/yyyy';
        topfig = [];  % image above animation axis
        UserData
    end
    properties (Dependent=true)
        hrange    % range of head contours
        prange    % range of stream function values for drawining stream lines
        crange    % value to be contoured for each of the UCN
        drange    % drawdown range
    end
    methods
        function o=animateObj(basename,varargin)
            %ANIMATEOBJ -- generates an animate object
            %
            % USAGE:
            %    animate = animateObj(basename,titleStrs,[,'head'][,'budget']})
            %    animate = animateObj(basename,{'tracer','temp'},'head','budget')
            %
            % titleStrs is a cell array holding the names of the species in
            % the simulation.
            % 'head'   indicates that heads should be loaded for simulation
            % 'budget' indicates that the budget file should be loaded for simulation
            % Concentrations, and heads can be simulated using methods of
            % animateObj, e.d.
            %  animate.concXS()
            %  animate.concYS()
            %  animate.concXY()
            %  animate.headXS()
            %  animate.headYS()
            %  animate.headXY()
            % TO 130615
                        
            if nargin==0
                return;
            end
            
            [o.t0   ,varargin] = getProp(varargin,{'t0','tstart'},0);
            [figPos,varargin ] = getProp(varargin,{'figpos','pos'},[]);
            [figName,varargin] = getProp(varargin,'fig','');

            if ~isempty(figName)
                error('%s: figName is not a property of %s',mfilename,mfilename);
            end                
            if ~isempty(figPos)
                error('%s: figPos is not a property of %s',mfilename,mfilename);
            end
            
            
            % basename is always necessary
            o.basename = basename;
            
            [o.gr,varargin] = getType(varargin,'gridObj',[]);

            % allow for keyworkd mask and patch to input a patch
            [o.patch ,varargin]  = getProp(varargin,'mask',o.patch);
            [o.patch ,varargin]  = getProp(varargin,'patch',o.patch);
            [o.psiMask,varargin] = getProp(varargin,'psiMask',o.psiMask);
            
            % check to see that patch as required fields
            if ~isempty(o.patch) && ...
                    (~isstruct(o.patch)     ||...
                     ~isfield( o.patch,'x') ||...
                     ~isfield( o.patch,'y') ||...
                     ~isfield( o.patch,'color'))
                error('%s: mask or patch must contain fields ''x'', ''y'' and ''color''',...
                    mfilename);
            end            
            % Is should not matter whether titleStrs is given as a single
            % cell array or as a list of strings
            k=0;
            [titleStrs,varargin] = getType(varargin,'cell',[]);
            if isempty(titleStrs)
                titleStrs = [varargin(:)'];
            end
            
            if isempty(titleStrs)
                error(['%s: animateObj constructor requires that you specify the basename\n',...
                      'followed by a set of titles of the parameters to be read and visualized.\n',...
                      'for example:\n',...
                      'animate = animateObj(basename,''temp'',''salinity'',''heads'',''budget'')\n',...
                      'SEE ALSO help animateObj'],mfilename);
            end
            
            i = strmatchi({'head','hds'},titleStrs);
            if i(1)
                if numel(i)>1
                    error('%s: head... is not unique in call. See labels %s',...
                        mfilename,sprintfs(' <<%s>>',titleStrs(i)));
                end
  
                cprintf('Keywords','%s: Reading HEADS from <<%s>>\n',mfilename,[basename '.HDS']);
                o.H=readDat([basename '.HDS']); o.H=maskHC(o.H,[-1000,1000],[NaN,NaN]);
                titleStrs(i)=[];
            end
            
            i = strmatchi({'drawd','ddn'},titleStrs);
            if i(1)
                if numel(i)>1
                    error('%s: drawd... is not unique in call. See labels %s',...
                        mfilename,sprintfs(' <<%s>>',titleStrs(i)));
                end
                cprintf('Keywords','%s: Reading DRAWDOWNS from <<%s>>\n',mfilename,[basename '.DDN']);
                o.D=readDat([basename '.DDN']); o.D=maskHC(o.D,[-1000,1000],[NaN,NaN]);
                titleStrs(i)=[];
            end
            
            i = strmatchi({'budg','bgt','balance'},titleStrs);
            if i
                if numel(i)>1
                    error('%s: budg... is not unique in call. See labels %s',...
                        mfilename,sprintfs(' <<%s>>',titleStrs(i)));
                end
                
                cprintf('Keywords','%s: Reading CELL BY CELL FLOWS from <<%s>>\n',mfilename,[basename,'.bgt']);
                
                %% Check if specific fields are to be read instead of all
                % format should be '...  'budget FLOWR FLOEWL' ...
                userLabels =  regexp(titleStrs{i},'\w+','match');
                if numel(userLabels)>1
                    o.B      = readBud([basename,'.BGT'],userLabels(2:end));  % get only flow rightface            
                else
                    o.B      = readBud([basename,'.BGT']);  % get only flow rightface            
                end
                o.B      = mf_Psi(o.B);
                o.B      = mf_setTime(o.basename,o.B);
                titleStrs(i)=[];
                
                o = o.maskPsi();
            end

            %% Remaining titles in titleStrs are the compoments from the
            %  transport model in order used by the transport model.
            %  This number may be less than or equal to the number of components
            %  actually used and in the order used MT3D00n.UNC where n is
            %  increasing and matched with the titleStrs.
            
            %%
            o.titleStr = titleStrs;      % names of components
            o.NCOMP = numel(titleStrs);  % number of components used
            
            fprintf('Components read: <<%s >>\n',sprintf(' ''%s''',titleStrs{:}));
                        
            % we know the names of the UNC files, so read them all
            for iComp=1:o.NCOMP
                try
                    cprintf('Keywords','%s: Reading <<%s>> from <<MT3D%03d.UCN>>\n',...
                        mfilename, upper(o.titleStr{iComp}),iComp);                    
                    o.UCN{iComp} = readMT3D(sprintf('MT3D%03d.UCN',iComp));
                catch ME
                    error(['%s: %s\n',...
                        'Can''t find or read file <<MT3D%03d.UCN>> for reading <<%s>>'],...
                        mfilename,ME.message,iComp,o.titleStr{iComp});
                end
                o.UCN{iComp} = maskHC(o.UCN{iComp},[0,Inf],[0,Inf]);
            end
            
            try
                o.time = [o.H.totim];
            catch ME
                fprintf('%s: %s\nH.totim not available, using time from concentration structs\n',...
                    mfilename,ME.message);
                try
                    o.time = [o.UCN{1}.time];
                catch ME
                    error('%s: %s\nNo compoments available to get time from.',...
                        mfilename,ME.message);
                end
            end
            
        end
        
        function prange = get.prange(o)
            if isempty(o.B)
                prange = [];
            else
                 prange = ContourRange(o.B,o.numberOfStreamLines,[],'Psi');
            end
        end
        function hrange = get.hrange(o)
            if isempty(o.H)
                hrange = [];
            else
                 hrange = ContourRange(o.H,o.numberOfHeadContours);
            end
        end
        function crange = get.crange(o)
            for iComp = o.NCOMP:-1:1
                if isempty(o.UCN{iComp})
                    crange{iComp} = [];
                else
                 crange{iComp} = [0 ContourRange(o.UCN{iComp},o.numberOfConcContours)];
                end
            end
        end
        function drange = get.drange(o)
            if isempty(o.D)
                drange = [];
            else
                 drange = ContourRange(o.D,o.numberOfHeadContours);
            end
        end
        function o = maskPsi(o,psiMask)
            %MASKSPI masks Psi on budget struct iwth given mask
            %
            % Example:
            %   o = o.maskPsi(psiMask)
            %   notice that size of psiMask must match size of o.B(it).Psi
            %
            % See also: maskHC maskBud
            %
            % TO 130504
            
            if nargin<2
                psiMask = o.psiMask;
            end
            if ~isempty(psiMask)
                if any(size(psiMask) ~= size(o.B(1).Psi))
                    error('%s: size of psiMask ([%d %d] must agree with size of Psi [%d %d]',...
                        mfilename,size(psiMask),size(o.psiMask));
                end
                o.psiMask = psiMask;
            end
            
            if isempty(o.psiMask)
                return;
            else
                for i=1:numel(o.B)
                    o.B(i).Psi(psiMask) = NaN;
                end
            end
        end
    end
end

        
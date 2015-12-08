classdef mnwObj
    %   TO 110804
    properties
        WELLID     % label SITE: textual well id
        TYPE       % WEL, WMN1, WMN2
        PRINT      % 0 or 1 whether or not to print this well in output files
        isaWEL     % treat as ordinary WEL instead of MNW
        Nr         % well number used to select Q and Conc columns in PER worksheet
        LOSSTYPE   % SKIN NONLINEAR ...
        NPER; NCOMP % mfLab number of stress periods and number of species components
        NX; NY; NZ % mfLab size of model (cells)
        x; ix      % mfLab x of well and position in grid
        y; iy      % mfLab y of well and position in grid
        z; iz      % mfLab z of well and position in grid
        zscrtop    % mfLab ztop   elevation(s)  of screen
        zscrbot    % mfLab bottom elevation(s)  of screen
        zGr        % mfLab layer elevation at well location
        zt; zb     % mfLab top and bottom of aquifer cells
        dz         % mfLab thickness of aquifer cells
        zm         % mfLab cell center elevation aquifer cells at well location
        LAYCBD     % mfLab whether or not confining bed below layers
        Qdes       % mfLab desired injection rate(if positive) of this well all stress periods
                   % specify only the first node of an MNW
        Conc       % mfLab injection conc this well all stress periods and all components
        rw         % MNW1 and MNW2 (MNW1: if rw<0 absolute value of cell to well conductance)
        rskin      % MNW2
        kskin      % MNW2
        skin       % MNW1, computed in setfQ, needs HK to be computed
        lossC      % MNW1 conefficient in non-linear cell to well constant
        lossP      % MNW1 power in non linear well loss function
        Hlim       % MNW1 useless parameter (limtis level in well, separtely for injection and extraction)
        Href       % MNW1 useells parameter but has to be set in MNW1 (set to Inf unless DD is on, then set
                   % to zero or so and set Hlim to Inf (Hopeless
                   % definitions in this manual)
        QWval      % MNW1 conc injected fluid. Specify for first not of MNW1 well only (needed in MODFLOW-GWT)
        DD         % MNW1 drawdown used or not ?
        Iqwgrp     % MNW1 water quality group identifyer, wells with same value are mixed unless QWval<0
        QCUT       % MNW1 flag (value QCUT: or Q-%CUT: isgnaling whether Qfrcmn and Qfrcmx are interpreted as
                   % rate or as a percentage of Qdes.
        Qfrcmn     % MNW1 min pumping rate to remain active
        Qfrcmx     % MNW1 min rate to reactivate the well
        
        %% Grid info
        LRC        % mfLab Layer Row Col of screened grid cells
        Idx        % mfLab Global index of screened grid cells
        L          % length of cell intersections
        fQ         % extracted fraction based on HK(Idx)*L/sum(HK(Idx)*L)
        P          % cell info
        %%         % user
        remark
        created
        UserData

    end
    properties (Constant=true)
        NNODES   = -1;    % only vertical MNW (multinode wells
        PUMPLOC  = 0;     % only default pomploc implemented
        QLIMIT   = 0;     % no limitations implemented
        PPFLAG   = 0;     % partial penetration flag
        PUMPCAP  = 0;     % no pumpcap constraints imiplemented
        ITYPEWEL = 2;     % itype WEL according to manuel
        ITYPEMNW = 27;    % itype MNW according to manuel
    end
    methods   
        function obj=mnwObj(wellid,type,xGr,yGr,zGr,LAYCBD,mnwnams,mnwvals,Qdes,Conc,HK)
            % wellnm is given explicitly to keeop mnwvals all numeric
            % zGr is either full 3D array of layer elevations or that of the
            % geven well position (or of uniform layer elevations)
            % zGr includes confining beds
            % LAYCBD is 1 for layers with a confining bed below them and zero without
            % mnwnams = header of columns of mnwvals which must include
            % Nr x y z rw rskin kskin
            % Qdes is injections for this MNW all stress periods Q(NPER,1);
            % Conc is injection concentration of this MNW all stress periods
            % and components Conc(NPER,NCOMP);
            
            if nargin==0 % genrate a single mnwObj 
                return;
            end
            
            if nargin==1  % generate Nx1 array of empty mnwObj
                N=wellid;
                obj(N)=mnwObj();
                return
            end
 
            obj.Nr     = mnwvals(strmatchi('Nr',mnwnams));

            obj.WELLID = wellid;  % textual id of well
            obj.TYPE   = type;    % textual type of well (WEL, WMN1, WMN2)
            obj.isaWEL = strcmp(obj.TYPE,'WEL');
 
            ltype = mnwvals(strmatchi('LOSSTYPE',mnwnams));
            if strcmp(obj.TYPE,'MNW1'),
                ltype = -abs(ltype);
            else
                ltype = +abs(ltype);
            end
            
            switch ltype
                case -3, obj.LOSSTYPE = 'NONLINEAR'; % MNW1
                case -2, obj.LOSSTYPE = 'LINEAR';    % MNW1
                case -1, obj.LOSSTYPE = 'SKIN';      % MNW1
                case  0, obj.LOSSTYPE = 'NONE';      % MNW2
                case  1, obj.LOSSTYPE = 'NONE';      % MNW2
                case  2, obj.LOSSTYPE = 'THIEM';     % MNW2
                case  3, obj.LOSSTYPE = 'SKIN';      % MNW2
                case  4, obj.LOSSTYPE = 'GENERAL';   % MNW2
                case  5, obj.LOSSTYPE = 'SPECIFYCWC';% MNW2
                otherwise
                    error(['Illegal LOSSTYPE for %s, use\n',...
                    '-1 = SKIN      (MNW1)\n', ...
                    '-2 = LINEAR    (MNW1)\n', ...
                    '-3 = NONLINEAR (MNW1)\n'...
                    ' 0 = NONE      (MNW1)\n'...
                    ' 1 = NONE      (MNW2)\n',...
                    ' 2 = THIEM     (MNW2)\n',...
                    ' 3 = SKIN      (MNW2)\n', ...
                    ' 4 = GENERAL   (MNW2)\n', ...
                    ' 5 = SPECIFYCWC(MNW2)\n',obj.TYPE]);
            end
            
            obj.x  = mnwvals(strmatchi('x',mnwnams)); % one or more values
            obj.y  = mnwvals(strmatchi('y',mnwnams)); % one or more values
            obj.z  = mnwvals(strmatchi('z',mnwnams)); % must be two values screen top and bot

            if rem(length(obj.z),2)
                error('%s: need an even number of z values to deduce tops and bottoms of screens in the same MNW (%d)',class(obj),length(obj.z));
            end
    
            % make sure this works for an arbitrary number of wells
            obj.zscrtop=obj.z(1:2:end); obj.zscrtop=obj.zscrtop(~isnan(obj.zscrtop));
            obj.zscrbot=obj.z(2:2:end); obj.zscrbot=obj.zscrbot(~isnan(obj.zscrbot));

            % extend x,y coordinates to length of screen so that it will
            % work for any number of vertical or non vertical screens, the
            % latter having bottom xy different from top xy 
            obj.x(end:length(obj.zscrtop))  =obj.x(end);
            obj.y(end:length(obj.zscrbot))  =obj.y(end);
            
            obj.rw = mnwvals(strmatchi('rw',mnwnams));
            [obj.ix,obj.iy,obj.iz] = xyzindex(obj.x,obj.y,obj.z,xGr,yGr,zGr);
            if any(isnan([obj.ix,obj.iy,obj.iz]))
                error('mnwObj:checkWells',...
                    'wells must be entirely within the model')
            end
            obj.LAYCBD = LAYCBD;
            obj.NX=length(xGr)-1;
            obj.NY=length(yGr)-1;
            obj.NZ=length(LAYCBD);

            if size(zGr,1)==1 && size(zGr,2)==1
                obj.zGr=zGr(1,1,:);
            else
                obj.zGr=zGr(obj.iy,obj.ix,:);
            end

            % Find the index, thickness and center of the cells
            % taking LAYCBD into consideration
            obj.dz=NaN(1,1,obj.NZ);
            obj.zm=NaN(1,1,obj.NZ);
            obj.zt=NaN(1,1,obj.NZ);
            obj.zb=NaN(1,1,obj.NZ);
            k=1;
            for izz=1:obj.NZ
                obj.zt(1,1,izz)=obj.zGr(1,1,k);
                obj.zb(1,1,izz)=obj.zGr(1,1,k+1);
                obj.dz(1,1,izz)= obj.zt(1,1,izz)-obj.zb(1,1,izz);
                obj.zm(1,1,izz)=(obj.zt(1,1,izz)+obj.zb(1,1,izz))/2;
                if LAYCBD(izz)>0,
                    k=k+2;
                else
                    k=k+1;
                end
            end
            
            % ensure obj.x,obj.y and obj.z have same lengths (necessary in
            % case of vertical well where x,y is specified only once.
            obj.x(end:length(obj.z))=obj.x(end);
            obj.y(end:length(obj.z))=obj.y(end);
            
            obj.P=linegridObj([obj.x(:),obj.y(:),obj.z(:)],xGr,yGr,zGr,LAYCBD);
            
            obj.LRC=[[obj.P.iz]' [obj.P.iy]' [obj.P.ix]'];
            obj.Idx= [obj.P.I];
            obj.ix = [obj.P.ix];
            obj.iy = [obj.P.iy];
            obj.iz = [obj.P.iz];
            obj.L  = [obj.P.L];
                
 %           obj.LRC=cellIndices([obj.Idx],[obj.NY,obj.NX,obj.NZ],'LRC');
            obj.Qdes   = Qdes;   obj.NPER  = size(Qdes,1);
            obj.Conc   = Conc;   obj.NCOMP = size(Conc,2);

            % MNW2 specific
            obj.rw    = mnwvals(strmatchi('rw',   mnwnams,'exact'));

            try
                obj.PRINT = mnwvals(strmatchi('PRINT',mnwnams));
                if obj.PRINT>0,
                    obj.PRINT=1;
                else
                    obj.PRINT=0;
                end
            catch ME
                obj.PRINT = 1;
                fprintf('%s %s %s not defined for well %s\n',...
                    class(obj),ME.identifier,'PRINT',obj.PRINT,obj.WELLID);
                fprintf('value set to %s\n','PRINT');
            end
            
            try
                obj.rskin = mnwvals(strmatchi('rskin',mnwnams));
            catch ME
                obj.rskin = obj.rw;
                fprintf('%s %s %s not defined for well %s\n',...
                    class(obj),ME.identifier,'rskin',obj.WELLID);
                fprintf('value set to %f\n',obj.rskin);
            end
            
            try
                obj.kskin = mnwvals(strmatchi('kskin',mnwnams));
            catch ME
                obj.kskin = 1e6; % make it huge
                fprintf('%s %s %s not defined for well %s\n',...
                    class(obj),ME.identifier,'kskin',obj.WELLID);
                fprintf('value set to %f\n',obj.kskin);
            end
            
            %% Specific for NMW1
            try
                obj.skin  = mnwvals(strmatchi('Skin',  mnwnams));
            catch ME,
                obj.skin = 0;  % set to no friction
                fprintf('%s %s %s not defined for wel %s\n',...
                class(obj),ME.identifier,'skin',obj.WELLID);
                fprintf('value set to %f\n',obj.skin);
            end
            
            try
                obj.QWval = mnwvals(strmatchi('QWval',  mnwnams));
            catch ME,
                obj.QWval=NaN;
                fprintf('%s %s %s not defined for well %s\n',...
                    class(obj),ME.identifier,'QWval',obj.WELLID);
                fprintf('value set to %f\n',obj.QWval);
            end
            
            try
                obj.Hlim  = mnwvals(strmatchi('Hlim',   mnwnams));
            catch ME,
                obj.Hlim =0;
                fprintf('%s %s %s not defined for well %s\n',...
                    class(obj),ME.identifier,'Hlim',obj.WELLID);
                fprintf('value set to %f\n',obj.Hlim);
            end
            
            try
                obj.Href  = mnwvals(strmatchi('Href',   mnwnams));
            catch ME,
                obj.Href =1e8; fprintf('setting Href in MNW1 to %g\n',obj.Href);
                fprintf('%s %s %s not defined for well %s\n',...
                    class(obj),ME.identifier,'Href',obj.WELLID);
                fprintf('value set to %f\n',obj.Href);
            end
            
            try
                obj.DD    = mnwvals(strmatchi('DD',     mnwnams,'exact'));
                switch obj.DD
                    case 0
                        obj.DD='';
                    case 1
                        obj.DD='DD';
                    otherwise
                        error('mnwObj:DD',...
                            '%s: DD must be wither 0 or 1, not %f\n',...
                            class(obj),obj.DD);
                end
            catch ME,
                obj.DD = 1;
                fprintf('%s %s %s not defined for well %s\n',...
                    class(obj),ME.identifier,'DD',obj.WELLID);
                fprintf('value set to %s\n',obj.DD);
            end
            
            try
                obj.lossC = mnwvals(strmatchi('lossC',mnwnams,'exact'));
            catch ME
                obj.lossC = 0;
                fprintf('%s %s %s not defined for well %s\n',...
                    class(obj),ME.identifier,'well loss factor lossC',obj.WELLID);
                fprintf('value set to %d\n',obj.lossC);
            end
            
            try
                obj.lossP = mnwvals(strmatchi('lossP',mnwnams,'exact'));
            catch ME
                obj.lossP = 0;
                fprintf('%s %s %s not defined for well %s\n',...
                    class(obj),ME.identifier,'well loss factor P',obj.WELLID);
                fprintf('value set to %d\n',obj.lossP);
            end
            
            try
                obj.Iqwgrp= mnwvals(strmatchi('Iqwgrp', mnwnams));
            catch ME,
                obj.Iqwgrp=obj.Nr;
                fprintf('%s %s %s not defined for well %s\n',...
                    class(obj),ME.identifier,'Iqwgrp',obj.WELLID);
                fprintf('value set to %d\n',obj.Iqwgrp);
            end
            
            try
                obj.QCUT = mnwvals(strmatchi('QCUT',mnwnams));
                if obj.QCUT>0, obj.QCUT=1; end
                switch obj.QCUT
                    case 0
                        obj.QCUT='Q-%CUT';
                    case 1
                        obj.QCUT='QCUT';
                    otherwise
                        throw(ME);
                end
            catch ME,
                obj.QCUT=NaN;
                fprintf('%s %s %s not defined for well %s\n',...
                class(obj),ME.identifier,'QCUT resp. Q-%CUT',obj.WELLID);
                fprintf('value set to %s\n',obj.QCUT);
            end
            
            try
                obj.Qfrcmn = mnwvals(strmatchi('Qfrcmn',mnwnams));
            catch ME
                obj.Qfrcmn=0;
                fprintf('%s %s %s not defined for well %s\n',...
                    class(obj),ME.identifier,'Qfrcmn',obj.WELLID);
                fprintf('value set to %d\n',obj.Qfrcmn);
            end
            try
                obj.Qfrcmx = mnwvals(strmatchi('Qfrcmx',mnwnams));
            catch ME
                obj.Qfrcmx=0;
                fprintf('%s %s %s not defined for well %s\n',...
                    class(obj),ME.identifier,'Qfrcmx',obj.WELLID);
                fprintf('value set to %d\n',obj.Qfrcmx);
            end

            %% Now that we have read in all possible data we may just as well
            %  Compute the skin also for ordinary wells (WEL) using kskin and
            %  rskin. This is done here according to Halford and Hanson
            %  (2002)
            
            obj.fQ=HK(obj.Idx(:)).* obj.L(:) / sum(HK(obj.Idx(:)).*obj.L(:));
            % because regular wells are recognized by their type being WEL
            % in the spreadsheet MNW, we compute it according to the MNW manual
            % (however, assuming isotropy and using only HK
            
            if strcmp(obj.TYPE,'WEL')
                obj.skin = ((HK(obj.Idx)./obj.kskin)-1)*log(obj.rskin/obj.rw);
            end
        end
        
        %% Writing MNW2 file
        function writeMNW(obj,fp)
            fprintf(fp,'%20s%10d  ; WELLID  NNODES\n',obj.NNODES);
            fprintf(fp,'%-20s%10d%10d%10d%10d  ; LOSSTYPE, PUMPLOC, Qlimit, PPFLAG, PUMPCAP\n',...
                obj.LOSSTYPE,obj.PUMPLOC,obj.Qlimit,obj.PPFLAG,obj.PUMPCAP);
            fprintf(fp,'%12f %12f %12f ;  rw rskin kskin\n',obj.rw,obj.rskin,obj.kskin);
        end
        function writeMNW_Q(obj,fp)
            for iper=1:length(obj.Qdes)
                fprintf(fp,'%20s%12f  ; WELLID  Qdes\n',obj.WELLID,obj.Q(iper));
            end
        end
        
        
        function obj=plot(obj,clr,lw,varargin)
            % USAGE mnwObj=plot(obj,clr,linewidth,vaargin)
            if nargin<3, lw=1; end
            if nargin<2, clr='b'; end
            for i=1:length(obj.zscrtop)
                plot3(obj.x,obj.y,[obj.zscrtop(i),obj.zscrbot(i)],[clr 'o-'],'linewidth',lw);
            end
            switch numel(varargin)
                case 0, text(obj.x(1),obj.y(1),obj.zscrtop(1),['  ' obj.WELLID]);
                case 2, text(obj.x(1),obj.y(1),obj.zscrtop(1),['  ' obj.WELLID],varargin{1},varargin{2}); 
                case 4, text(obj.x(1),obj.y(1),obj.zscrtop(1),['  ' obj.WELLID],varargin{1},varargin{2},varargin{3},varargin{4});
                case 6, text(obj.x(1),obj.y(1),obj.zscrtop(1),['  ' obj.WELLID],varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6});
                case 8, text(obj.x(1),obj.y(1),obj.zscrtop(1),['  ' obj.WELLID],varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6},varargin{7},varargin{8});
                otherwise
                    error('too many or odd number of input arguments in mnwObj.plot beyond 3');
            end
        end
        function obj=plot2d(obj,clr,lw,varargin)
            % USAGE mnwObj=plot(obj,clr,linewidth,vaargin)
            if nargin<3, lw=1; end
            if nargin<2, clr='b'; end
            for i=1:length(obj.zscrtop)
                plot([obj.P.x],[obj.P.z],[clr '-'],'linewidth',lw);
            end
            switch numel(varargin)
                case 0, text(obj.x(1),obj.zscrtop(1),['  ' obj.WELLID]);
                case 2, text(obj.x(1),obj.zscrtop(1),['  ' obj.WELLID],varargin{1},varargin{2}); 
                case 4, text(obj.x(1),obj.zscrtop(1),['  ' obj.WELLID],varargin{1},varargin{2},varargin{3},varargin{4});
                case 6, text(obj.x(1),obj.zscrtop(1),['  ' obj.WELLID],varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6});
                case 8, text(obj.x(1),obj.zscrtop(1),['  ' obj.WELLID],varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6},varargin{7},varargin{8});
                otherwise
                    error('too many or odd number of input arguments in mnwObj.plot beyond 3');
            end
        end
        function WEL=wel(obj)
            NScrCell=length(obj.Idx);
            u=ones(NScrCell,1);
            WEL=NaN(obj.NPER*NScrCell,5);
            for iper=1:obj.NPER
                WEL((iper-1)*NScrCell+(1:NScrCell),:)=...
                    [iper*u obj.LRC obj.fQ*obj.Qdes(iper,1)];
            end
        end
        function PNTSRC=pntsrc(obj)
            if ~obj.isaWEL,          % it knows what to print for: WEL or MNW
                ITYPE=obj.ITYPEMNW;
%                NScrCell=length(obj.Idx(1));  % < only one per multinode well
                NScrCell=length(obj.Idx);  % < only one per multinode well
            else
                ITYPE=obj.ITYPEWEL;
                NScrCell=length(obj.Idx);     % < take all for ordinary bdl
            end
            PNTSRC=NaN(obj.NPER*NScrCell,6+obj.NCOMP);
            u=ones(NScrCell,1);
            for iper=1:obj.NPER
                PNTSRC((iper-1)*NScrCell+(1:NScrCell),:)=...
                    [iper*u obj.LRC(1:NScrCell,:) u*obj.Conc(iper,1) u*ITYPE u*obj.Conc(iper,:)];
            end
        end
        function obj=set.TYPE(obj,wellid)
            switch wellid
                case {'WEL', 'WELL'},      obj.TYPE='WEL';
                case {'MNW1', 'MN1' 'M1'}, obj.TYPE='MNW1';
                case {'MNW2', 'MN2','M2'}, obj.TYPE='MNW2';
                otherwise
                    error('%s TYPE must be one of WEL MNW1 MNW2\n',class(obj));
            end
        end
    end
end


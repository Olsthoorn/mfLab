classdef piezObj
    properties
        test
        aspect
        include
        name
        r
        z
        DDN
        use
        UserData
    end
    methods
        function o = piezObj(basename,piezSheet,includeSheet,dataSheet)
            %%PIEZOBJ -- generates a piezometer object
            %
            % USAGE:
            %   piez = piezObj(basename,piezSheetNm,includeSheetNm,dataSheetNm)
            %
            %   piezSheetNm contains table with fields [name	r	z]
            %   includeSheetNm contains table with fields [test	name
            %   include]
            %   dataSheetNm contains table with fields
            %          [test	name	piezom	rdum	tmin	ddn]
                        
            if nargin==0
                return;
            end
            
            [pzHdr ,pzVals ,pzTxtHdr , pzTxtVals] = getExcelData(basename,piezSheet,'Hor');
            pz_nmCol     = strmatchi('name'   ,pzTxtHdr);
            pz_zCol      = strmatchi('zm'     ,pzHdr,    'exact');
            pz_rCol      = strmatchi('r'      ,pzHdr,    'exact');
            
            if nargin>2
                [inclHdr,inclVals,inclTxtHdr,inclTxtVals] = getExcelData(basename,includeSheet,'Hor');
                incl_nmCol   = strmatchi('name'   ,inclTxtHdr);
                incl_inclCol = strmatchi('include',inclHdr);
                incl_tstCol  = strmatchi('test'   ,inclHdr);
            end            
            
            if nargin>3
                [ddnHdr,ddnVals,ddnTxtHdr,ddnTxtVals] = getExcelData(basename,dataSheet,'Hor');
                ddn_nmCol    = strmatchi('name'   ,ddnTxtHdr');
                ddn_tstCol   = strmatchi('test'   ,ddnHdr');
                ddn_ddnCol  = strmatchi({'tmin','ddn'},ddnHdr);
            end

            for i = size(inclVals,1):-1:1
                o(i).name    = inclTxtVals(i,incl_nmCol);
                o(i).name    = o(i).name{1};
                
                if nargin>2
                    o(i).test    = inclVals(i,incl_tstCol);
                    o(i).include = inclVals(i,incl_inclCol);
                end
                j = strmatchi(o(i).name,pzTxtVals(:,pz_nmCol),'exact');

                o(i).r    =  pzVals(j,pz_rCol);
                o(i).z    =  pzVals(j,pz_zCol);

                if nargin>3
                    I = strcmp(o(i).name,ddnTxtVals(:,ddn_nmCol)) & ...
                        ddnVals(:,ddn_tstCol)==o(i).test;
                    o(i).DDN = ddnVals(I,ddn_ddnCol);
                    o(i).DDN(:,1) = o(i).DDN(:,1)/(60*24);
                    o(i).use = true(size(o(i).DDN(:,1)));
                end
            end
        end
        function o = plot(o,varargin)
            [ax,varargin] = getProp(varargin,'axis',[]);
            if isempty(ax)
                hold on;
                [ax,varargin] = getNext(varargin,'axis',gca);
            end
            [testNrs,varargin] = getProp(varargin,'test',[]);
            if ~isempty(testNrs)
                I = ismember([o.test],testNrs) & [o.include];
            else
                [piezNms,~] = getProp(varargin,'piez',{});
                if ~isempty(piezNms)
                    if ~iscell(piezNms), piezNms={piezNms}; end
                    I = ismember({o.name},piezNms) & [o.include];
                else
                    I = [o.include]~=0;
                end
            end
            I = find(I);
            for i=numel(I):-1:1
                if ~any(o(I(i)).use)
                    I(i)=[];
                end
            end
            if isempty(I)
                fprintf('Nothing to plot !\n');
                delete(gcf);
                return;
            end
            leg={};
            leg{1,numel(I)}='';
            hold on
            k=0;
            for i=1:numel(I)
                j = I(i);
                k=k+1;
                leg{k} = sprintf('%s z=%.2f r=%.1f',o(j).name,o(j).r,o(j).z); %#ok
                plot(ax,o(j).DDN(o(j).use,1),o(j).DDN(o(j).use,2),['o' mf_color(j)]);
                if size(o(j).DDN,2)>2
                    k=k+1;
                    leg{k} = sprintf('%s z=%.2f r=%.1f, model',o(j).name,o(j).r,o(j).z); %#ok
                    plot(ax,o(j).DDN(o(j).use,1),o(j).DDN(o(j).use,3),mf_color(j));
                end
            end
            legend(ax,leg{:},2);
        end
        function o = useAspect(o,aspects)
            %%USEASPECT -- uses USE for using aspect
            % USAGE:
            %   piez = piez.useAspect([1 2 ...]); % set use for these aspects
            %   piez = piez.useAspect();          % set use to true for all times            
            % This resets the use to all true
            if nargin<2
                for io=1:numel(o)
                    o(io).use = true(size(o.DDN(:,1)));
                end
                return;
            end
            
            % using aspects implies use == false for all not use of it
            for io = 1:numel(o)
                % always start setting use to false
                o(io).use = false(size(o(io).DDN(:,1)));
                
                tmin = o(io).DDN(end,1);
                tmax = o(io).DDN(   1,1);
                for ia = 1:numel(aspects)
                    for j = 1:numel(o(io).aspect)
                        if o(io).aspect(j).nr == aspects(ia)
                            tmin = min(tmin,o(io).aspect(j).timeSpan(  1));
                            tmax = max(tmax,o(io).aspect(j).timeSpan(end));
                        end
                    end
                end
                o(io).use = o(io).DDN(:,1)>=tmin & o(io).DDN(:,1)<=tmax;
            end
        end
        function o = setAspect(o,basename,sheetName)
            %%SETASPECT -- sets a calibration aspect defined by aspect
            % number and time span. One can include any or all aspects or
            % none. None = all data. Including aspecs implies using the
            % only the data within the time span belonging to the aspect.
            % This facilitates calibrating certain aspects of the problem.
            [aspHdr,aspVals,aspTxtHdr,aspTxt] = getExcelData(basename,sheetName,'Hor');
            testCol = strmatchi('test',  aspHdr);
            nameCol = strmatchi('name',  aspTxtHdr);
            remCol  = strmatchi('remark',aspTxtHdr);
            nrCol   = strmatchi('aspect',aspHdr);
            tCol    = strmatchi('time', aspHdr);  % yields start and end time 
            for io=1:numel(o)
                o(io).use = true(size(o(io).DDN(:,1)));
            end
            for itest = 1:max([o.test])
                I = find([o.test] == itest);
                A = aspVals(aspVals(:,testCol)==itest,:);
                T = aspTxt( aspVals(:,testCol)==itest,:);
                for ip = numel(I):-1:1
                    J = find(ismember(T(:,nameCol),o(I(ip)).name));
                    for j=numel(J):-1:1
                       o(I(ip)).aspect(j).nr = A(J(j),nrCol);
                       o(I(ip)).aspect(j).timeSpan = A(J(j),tCol);
                       o(I(ip)).aspect(j).remark   = T{J(j),remCol};
                    end
                end
            end
        end
        function o = hantMod(o,Q,k,Ss,screen,D)
            %%HANTUSHEMOD -- compute modified Hantush solution
            %
            % USAGE piez = piez.hantushEmod(hk,kv,Ss,well)
            %     hk = hor k
            %     vk = vert k
            %     Ss = specific storage coefficient
            %     well= single well object contaning z of screen and Q
            
            for i=numel(o):-1:1
                t = o(i).DDN(:,1);
                o(i).UserData.theis        = hantushMod(Q,k,Ss,o(i).r,t,D);
                o(i).UserData.hantModEarly = hantushMod(Q,k,Ss,o(i).r,o(i).z,t,screen);
                o(i).UserData.hantModLate  = hantushMod(Q,k,Ss,o(i).r,o(i).z,t,screen,D);
                o(i).UserData.brug         = hantushMod(Q,k,Ss,o(i).r,o(i).z,t,mean(screen),D);
                o(i).UserData.brugInf      = hantushMod(Q,k,Ss,o(i).r,o(i).z,t,mean(screen));
            end
        end
        
    end
end
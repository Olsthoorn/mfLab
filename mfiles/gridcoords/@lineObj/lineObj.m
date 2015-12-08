classdef lineObj < pointObj
%LINEOBJ: generate line objects, defined in a spreadsheet
% lineObj objects are used to define line boundary conditions like rivers,
% line drains, line GHB, line CHD boundaries and lines with specified head
% and flux.
% lineObj is a subclass of poinObj and a superclass of area2Obj
%
% USAGE:
%       var  = lineObj([gr,] basename,sheetNm,type[,options[,plane[,plane[,plane]]]])
%
% EXAMPLES:
%       any = lineObj([gr,]basename,lineObjSheetNm,'NIL')
%       riv = lineObj([gr,]basename,lineObjSheetNm,'RIV',DEM,[],DEM)
%       drn = lineObj([gr,]basename,lineObjSheetNm,'DRN',DEM}
%       ghb = lineObj([gr,]basename,lineObjSheetNm,'GHB',DEM}
%       chd = lineObj([gr,]basename,lineObjSheetNm,'CHD',DEM,DEM)
%       flux= lineObj([gr,]basename,lineObjSheetNm,'WEL',MultiplierArrayOfSizeNyNx)
%
% see extended explantion with pointObj (help pointObj)
%
% Notice, yo can put the gr anywhere in the argument list
%
% SEE ALSO wellObj pointObj lineObj area2Obj gridObj.line
%
% TO 130818
    
    properties
        cellLineVals  % same as vertex but for the cells [P_idx]: [NcellOnPolygon ,size(vertex,o(io).iColx)]
        P             % gridLineObjects, see gridLObj (check)
        L             % cells intersected by each line (cell array)
        sLine         % line length along vertices
        sCell         % line length to mid of intersection along intersected cells

    end
    properties (Dependent=true)
        length        % length of the lineObject
    end
    methods
        function o=lineObj(varargin)
            %LINEOBJ: constructor of lineObj
            %
            % TO 130818
            
            o = o@pointObj(varargin{:});
            
            if nargin<2, return; end
            
            isLineObj = isa(o,'lineObj'); % else area2Obj
 
            %% Get workbookName, sheetName and bcn type
            
            %% Process input in the same manner as in super class to make sure
            %  we can reprocess some of the argument in a subclass-specific way
            [~,varargin] = getProp(varargin,'active','');
            [~,varargin] = getProp(varargin,'conc',[]);
            [~,varargin] = getProp(varargin,'struct',{});
            
            [gr,varargin] = getType(varargin,'gridObj',[]);
            [~ ,varargin] = getNext(varargin,'char',[]);  % basename
            [~ ,varargin] = getNext(varargin,'char',[]);  % sheetname
            [~ ,varargin] = getNext(varargin,'char',[]);  % type

            [colHdr_  ,varargin] = getNext(varargin,'char',''); %colHdr
            if ~isempty(colHdr_)
                [~ ,varargin] = getNext(varargin,{'char','double','cell'},[]);
            end
            
            if isempty(gr), return; end
            
            %% ======== Put lines into grid =====================
            msgId = 'lineObj:tooFewVertices';

            for io = numel(o):-1:1
                
                if strcmpi(o(io).name,'north')
                    fprintf \n
                end
                
                jx   = o(io).iColx;
                jy   = o(io).iColy;
                
                %% Intersect the lineObj with the model
                % For that we need the elevations of te vertices of the
                % line. There is a precedence: absolute z, relative z or
                % ground surface:
                
                % already in pointObj file we made sure every object has
                % its z or zRel column
                if size(o(io).vertex,1) < 2
                    warning('on','m:m');
                    warning('m:m','%s <<%s>> has only one vertex, will be removed',class(o(io)),o(io).name);
                    warning('off','m:m')
                    o(io)=[];
                    continue;  % silently remove
                end

                switch o(io).type
                    case {'NIL','WEL','FLUX','CHD','GHB'}
                        jz   = strmatchi('z'   ,o(io).vertexHdr,'exact');
                        if jz % use 3D z-coord for elevation of object                                        
                            [o(io).P,o(io).L] = gr.lineObjects([o(io).vertex(:,[jx jy]) o(io).vertex(:,jz)]);
                        else  % use zrel column for elevation                                   
                            jzR  = strmatchi('zRel',o(io).vertexHdr); 
                            if jzR
                                [o(io).P,o(io).L] = gr.lineObjects( o(io).vertex(:,[jx jy]),o(io).vertex(:,jzR));
                            end
                        end
                    case {'DRN','RIV'}
                        % use preferred column given order in which the column
                        % names are specified in prefs below. We'll take
                        % the first one that is present in the hdrs
                        % to determine the z-coordinates of the object
                        hdrs  = lower(o(io).vertexHdr);
                        if strcmpi('RIV',o(io).type)
                            prefs = lower({'z','zB','zRel','Hb','Db'});
                        else
                            prefs = lower({'z','zB','zRel','Hb','Db','H','D'});
                        end
                        i  = find(ismember(prefs,hdrs),1,'first');
                        jz = find(ismember(hdrs,prefs{i}),1,'first');
                        switch i
                            case {1,2,4,6} % absolute z given
                               [o(io).P,o(io).L] = gr.lineObjects([o(io).vertex(:,[jx jy]) o(io).vertex(:,jz)]);
                            case {5,7} % Db or D given
                                z = interp2(gr.Xc,gr.Yc,gr.Z(:,:,1), o(io).vertex(:,jx), o(io).vertex(:,jy))  - o(io).vertex(:,jz);
                                [o(io).P,o(io).L] = gr.lineObjects( [o(io).vertex(:,[jx jy]),z]);
                            case 3  %zRel given
                                [o(io).P,o(io).L] = gr.lineObjects( o(io).vertex(:,[jx jy]),o(io).vertex(:,jz));
                            otherwise
                        % skip                            
                        end
                end
                
                if numel(o(io).P)<2
                    warning('on',msgId);
                    warning(msgId,['%s %s has %d points (<2) will be removed !!.\n',...
                        'Check that its points are within the model boundaries'],...
                        class(o(io)),o(io).name,numel(o(io)));
                    o(io)=[];
                    warning('off',msgId);
                    continue;
                end

                
                o(io).Idx = [o(io).P.idx];
                o(io).A   = [o(io).P.L];    % to compute specific discharge
                
                % Notice: properties of P are
                % x y z zR     ix iy iz   idx  L  xm ym zm

                %% Interpolate all values to their grid locations
 
                % cumulative line length at vertices
                line = o(io).vertex(:,[jx jy]);
                o(io).sLine = [0; cumsum(sqrt(sum(diff(line,1,1).^2,2)))];

                %Cumulative length along the intersected cells
                % to the mid of the intersection with the cell
                line = [[o(io).P.xm]; [o(io).P.ym]].';
                o(io).sCell = [0; cumsum(sqrt(sum(diff(line,1,1).^2,2)))];
                % elevation is known to P don't need to compute it

                % Interpolate the remaining values in the vertices to the cell centers
                % make sure the order of the values remains as they are in the
                % worksheet
                o(io).cellLineVals = interp1(o(io).sLine,o(io).vertex,o(io).sCell);
                
                %% processing further command-line arguments.
                % These are assumed to be planes of size [Ny,Nx]
                
                % Set the planes and select the cells from them
                % depending on whether we have a lineObj or an area3Obj
                
                if isLineObj
                    for iPlane = 1:min(o(io).narg,numel(varargin))
                        if ~isempty(varargin{iPlane})
                            % if specified it must be a plane of the size of a layer
                            if ~all(size(varargin{iPlane}(:,:,1))==[gr.Ny,gr.Nx])
                                error('%s: argument must be of size Ny,Nx = [%d,%d]',mfilename,gr.Ny,gr.Nx);
                            end

                            % select the cells from the plane
                            % ensure we index in a 2D plane only:
                            Nxy = prod(o(io).grSize(1:2)); % cells in plane
                            Idx = rem([o(io).P.idx]-1,Nxy)+1;
                            o(io).V{iPlane} = varargin{iPlane}(Idx);

                            % Guarantee column vector
                            o(io).V{iPlane} = o(io).V{iPlane}(:);                    
                        end
                    end
                end
            end
        end
        
        function [T,o] = kD(o,varargin)
            %LINOBJ/KD -- computes total transmissivity arcross the line
            %object, in fact a conductance that multiplied by the gradient
            %across this line yields the total flow.
            % The value is integrate along the profiel but separated per
            % layer.
            % USAGE: [T,o] = lineObj.kD(gr,HK,H);
            %
            %  T is transmissivity vector (Nlay x Nobj)
            %  o is the objects back with kD added to their UserData
            %  gr = gridObj
            %  HK = 3D horizontal k array
            %  H  = head struct as read with readDat, containgin the heads.
            %       H(end).values will be used for the comptutation of kD.
            % TO 131011
            
            [gr,varargin] = getType(varargin,'gridObj',[]);
            [H ,varargin] = getType(varargin,'struct',[]);
            [HK,~       ] = getType(varargin,'double',[]);
            if isempty(gr)
                error('Requires gridObj');
            end
            if isempty(H)
                error('Requires head struct as read by readDat');
            end
            if isempty(HK)
                error('Requries HK  (horizontal conductivity array');
            end
            
            T = NaN(gr.Nlay,numel(o));
            
            for io=numel(o):-1:1
                L_    = ones(gr.Nlay,1) * [o(io).P.L];
                k     = NaN(size(L_));
                D     = NaN(size(L_));
                zB    = NaN(size(L_));
                h     = NaN(size(L_));
                
                IdTop = gr.IdTop([o(io).P.idx]);
                for iLay=numel(gr.Nlay:-1:1)
                    idx = IdTop + gr.Nxy*(iLay-1);
                    k( iLay,:) = HK(   idx);
                    D( iLay,:) = gr.DZ(idx);
                    zB(iLay,:) = gr.Z(idx+gr.Nxy);
                    h( iLay,:) = H(end).values(idx);
                end

                D = min(max(0,h-zB),D);

                tr = k .* D .* L_; tr(isnan(tr)) = 0;
                T(:,io) = sum(tr,2);
                o(io).UserData.kD = T(:,io);
            end            
        end
        
        function fill2(o,varargin)
            %FILL2: fill a cross section
            %
            % USAGE lineObj.fill2([axis],gr[,clrList,plotOptions]);
            %    axis is handle to axis, default is gca
            %    gr is the gridObj of the model
            %    list of matlab colors like 'rgbmckyw'
            %    plotOptions such as 'lineWidth',3 in name,value pairs
            %
            % TO 130823
            clrs = 'brgkmcyw';
            [gr,  varargin] = getType(varargin,'gridObj',[]);
            [ax,  varargin] = getNext(varargin,'axis',gca);
            [clrList, ~   ] = getNext(varargin,'char',[]);
            if isempty(clrList) || ~isColor(clrList)
                clrList = clrs;
            else
                varargin(1)=[];
            end
            
            if isempty(gr)
                error('%s: a grid is required as input argument',mfilename);
            end
            
            for io=1:numel(o)
                % plot a profile

                % vals must be of size Ny,Nx,Nlay+1
                if any(gr.size ~= o(io).grSize)
                    error('%s: size vals [%d %d] ~= grSize [%d %d]',...
                        mfilename,gr.size,o(io).grSize(1:2));
                end

                % the field has to be of size gr.Z (elevations, not cells)
                s    = o(io).sCell(:)';

                xm=[o(io).P.xm];
                ym=[o(io).P.ym];

                for iL=1:gr.Nlay
                    zt = interp2(gr.xc,gr.yc,gr.Z(:,:,iL  ),xm,ym);
                    zb = interp2(gr.xc,gr.yc,gr.Z(:,:,iL+1),xm,ym);
                    fill([s s(end:-1:1)],[zt zb(end:-1:1)] ,mf_color(iL,clrList),varargin{:},'parent',ax);
                end
            end
        end
        
        function length = get.length(o)
            length = sum(sqrt(diff(o.vertex(:,o.iColx)).^2+diff(o.vertex(:,o.iColy)).^2));
        end

        function fill3(o,varargin)
            %FILL3: fill a cross section in 3D space
            %
            % USAGE lineObj.fill3([axis],gr[,clrList,plotOptions]);
            %    axis is handle to axis, default is gca
            %    clrList is a list of Matlab colors like 'brgkmcyw'
            %    plotOptions such as 'lineWidth',3 in name,value pairs
            %
            % TO 130823
            clrs = 'brgkmcyw'; % default colors
            
            [gr,  varargin] = getType(varargin,'gridObj',[]);
            [ax,  varargin] = getNext(varargin,'axis',gca);
            [clrList, ~       ] = getNext(varargin,'char',[]);

            if isempty(gr)
                error('%s: gridObj is required as input',mfilename);
            end
            
            if ~isempty(clrList) && isColor(clrList), 
                varargin(1)=[];
            else
                clrList = clrs;
            end
                   
            if isempty(varargin)
                varargin = {'lineWidth',1};
            end                

            for io=1:numel(o)
                x   = [o(io).P.xm];
                y   = [o(io).P.ym];
                for iL=1:gr.Nlay
                    z1  = interp2(gr.xc,gr.yc,gr.Z(:,:,iL)  ,x,y);
                    z2  = interp2(gr.xc,gr.yc,gr.Z(:,:,iL+1),x,y);
                    patch([x x(end:-1:1)],[y y(end:-1:1)], [z1 z2(end:-1:1)] ,mf_color(iL,clrList),varargin{:},'parent',ax);
                end
            end
            view(3);
        end

        function Qavg = mean(o)
            %%LINEOBJ/MEANQ --- computes mean Q voor each object over all times
            % requires that setQ method was run.
            if isempt(o(1).Q)
                o=o.setQ;
            end
            Qavg = NaN(size(o));
            for io=1:numel(o)
                Qavg = mean(o(io).Q);
            end
        end
        
        function hdl = plot(o,varargin)
            %PLOT: plots contour of line or area object to see its position
            %
            % USAGE lineObj.plot([axis],lineSpec[,plotOptions][,'crit',crit,'clrs',clrs[,B]]);
            %    axis is handle to axis, default is gca
            %    lineSpec is a legal line specifiiation combining a color and a line type
            %    like 'r-', 'b--'
            %    plotOptions such as 'lineWidth',3 in name,value pairs.
            %    'crit',crit   define criteria for the total net flow to
            %    decide what color to plot the lineObj in.
            %    'clrs',clrs defines the colors with clrs like 'rbgkmcy'.
            %    The criteria are compared against the total flow of the
            %    line objects, which must have been computed using setQ. If
            %    not, it can be computed with the current method by adding
            %    the buget struct B (as read by readBud) as one of the
            %    argument of the call.
            %
            % SEE ALSO: lineObj/plot2 lineObj/plot3 lineObj/fill
            %
            % TO 130823 131221
            
            [crit,varargin] = getProp(varargin,'crit',[]);
            [clrs,varargin] = getProp(varargin,{'color','clr'},'');
            [B   ,varargin] = getType(varargin,'struct',[]);
            [ax  ,varargin] = getType(varargin,'axis',gca);

             if ~isempty(crit)
                 crit = unique(crit);
                 if isempty(clrs)
                     clrs = 'rbgkmcyrbgkmcy'; % default colors;
                 else
                     dcrClr = numel(clrs)-numel(crit);
                     if  dcrClr == 0
                         clrs= [clrs 'k'];
                     elseif dcrClr>=1
                         clrs = clrs(1:numel(crit)+1);
                     else
                         error('you shlould specify the number of colors equal to the number of criteria + 2');
                     end
                 end
                 if ~isempty(B)
                     o = o.setQ(B);
                 elseif isempty(o(1).UserData.(['Q',o(1).type]))
                     error(['Can''t decide on colors because flow Q of this object was not set.\n',...
                         'Add the B read by readBud(..) as one of the arguments of the call.']);
                 end
             end
            
            hdl=[];
            for io=1:numel(o)
                if isempty(crit)
                    if isempty(varargin)
                        args = {[mf_color(io), '-']};
                    elseif ~isLineSpec(varargin{1})
                        args = ['b' varargin];
                    else
                        args = varargin;
                    end
                else
                    % if criterium is used for coloring object then corresponding color
                    % is used for objects with abs(Q) above crit
                    % example crit [-1000 0 1000],'clrs' 'rbg' 
                    if ~isempty(o(io).UserData.(['Q',o(io).type])) 
                        % color according to mean value of Q
                        Qmean = mean(sum(o(io).UserData.(['Q',o(io).type])));
                        iClr = ceil(interp1(unique([-1e20 crit 1e20]),0:numel(crit)+1,Qmean));
                        if isnan(iClr)
                            iClr=1;
                        end
                        args = [clrs(iClr) varargin];
                    end                
                end
                xm = [o(io).P.xm];
                ym = [o(io).P.ym];
                hdl(io) = plot(ax,xm,ym,args{:}); %#ok

            end
        end        
                  
        function plot2(o,varargin)
            %PLOT2: plots values along lineObj in 2D vertical section
            %
            % USAGE lineObj.plot([axis],field[,lineSpec,plotOptions]);
            %    axis is handle to axis, default is gca
            %    field is array of the size of a grid layer, usually DEM or
            %    heads like H(end).values.
            %    lineSpec is character of the lines (color, type)
            %    plotOptions such as 'lineWidth',3 in name,value pairs
            %
            % TO 130823
            
            [gr,  varargin] = getType(varargin,'gridObj',[]);
            [ax,  varargin] = getNext(varargin,'axis',gca);
            [vals,varargin] = getNext(varargin,'double',[]);
            for io=1:numel(o)
                if isempty(vals)
                    if ~isempty(varargin)
                        plot(ax,[o(io).sCell],[o(io).P.zm],varargin{:});
                    else
                        plot(ax,[o(io).sCell],[o(io).P.zm],o(io).EdgeColor);
                    end
                else
                    if any(size(vals(:,:,1)) ~=o(io).grSize(1:2))
                        error('%s: size vals [%d %d] ~= grSize [%d %d]',...
                            mfilename,sv,o(io).grSize(1:2));
                    end

                    if isempty(varargin)
                        varargin = {'lineWidth',1};
                    end
                                        
                    if exist('gr','var')
                        Iz=[o(io).P.iz];
                        xm=[o(io).P.xm];
                        ym=[o(io).P.ym];
                        z = NaN(size(Iz)); 
                        for iL=1:gr.Nlay
                            I = Iz==iL;
                            z(I) = interp2(gr.xc,gr.yc,vals(:,:,iL),xm(I),ym(I));
                        end
                        plot(ax,o(io).sCell,z,varargin{:});
                    else
                        plot(ax,o(io).sCell,vals([o(io).P.idx]),varargin{:});
                    end
                end
            end
        end

        function plot3(o,varargin)
            % PLOT3: plots values along lineObj
            %
            % USAGE lineObj.plot3([axis],field[,lineSpec,plotOptions]);
            %    axis is handle to axis, default is gca
            %    field is array of the size of a grid layer, usually DEM or
            %    heads like H(end).values.
            %    lineSpec is character of the lines (color, type)
            %    plotOptions such as 'lineWidth',3 in name,value pairs
            %
            % TO 130823
            
            [gr,  varargin] = getType(varargin,'gridObj',[]);
            [ax,  varargin] = getNext(varargin,'axis',gca);
            [vals,varargin] = getNext(varargin,'double',[]);
            for io=1:numel(o)
                if isempty(o(io).P)
                    fprintf(2,'WARNING: %s/plot3(), %s ''%s'', has empty P field. Remedy add gridObj to argument list.\n',...
                        mfilename,class(o(io)),o(io).name);
                    continue;
                end
                xm = [o(io).P.xm];
                ym = [o(io).P.ym];
                zm = [o(io).P.zm];

                % close area2Obj in this plot
                if isa(o(io),'area2Obj')
                    xm = [xm xm(1)]; %#ok
                    ym = [ym ym(1)]; %#ok
                    zm = [zm zm(1)]; %#ok
                end

                if isempty(varargin)
                    varargin{1} = o(io).EdgeColor;
                elseif ~isLineSpec(varargin{1})
                    varargin = [mf_color(io) varargin]; %#ok
                end

                if isempty(vals)
                    plot3(xm,ym,zm,varargin{:});
                    continue;
                end
                
                % check size of vals
                sv = size(vals);
                if length(sv)<3, sv = [sv 1]; end %#ok

                if any(sv ~=o(io).grSize)
                    error('%s: size vals [%d %d] ~= grSize [%d %d]',...
                        mfilename,sv,o(io).grSize(1:2));
                end
                
                Idx= [o(io).P.idx];
                if isa(o(io),'areaObj')
                    Idx = [Idx Idx(1)]; %#ok
                end
                
                if isempty(gr)
                    plot3(gr.xm,gr.ym,gr.zm,vals(Idx),varargin{:});
                    continue;
                end
                
                IL = floor((Idx-1)/gr.Nxy); % layer number - 1
                
                vm = NaN(size(o(io).P));
                for iLay=1:gr.Nlay
                    I = IL==iLay-1;
                    if any(I)
                        vm(I) = interp2(gr.xc,gr.yc,vals(:,:,iLay),xm(I),ym(I));
                    end
                end
                    
                if isa(o(io),'area2Obj') % close the polygon
                    vm = [vm; vm(1)]; %#ok
                end

                plot3(ax,xm,ym,vm,varargin{:});
            end
        end

        function o = check(o,IBOUND)
            %CHECK: remove cells that conincide with IBOUND==0
            mesId = [class(o) ':pointsAtZeroIBOUND'];
            warning('on',mesId);
            for io=numel(o):-1:1
                if nargin<2 || ~all( size(IBOUND(:,:,1))==o(io).grSize(1:2) ) || size(IBOUND,3)~=o(io).grSize(3)
                    error('%s: need IBOUND in call and its size must be %s',...
                        mfilename,sprintf(' %d',o(io).grSize));
                end
                nEl1 = numel(o(io).Idx);
                J = find(IBOUND(o(io).Idx)~=0);
                o(io).cellLineVals    = o(io).cellLineVals(J,:);
                o(io).P               = o(io).P(J);
                o(io).A               = o(io).A(J);
                o(io).sCell           = o(io).sCell(J);
                for i=1:numel(o(io).V)
                    if numel(o(io).V{i})>1
                        o(io).V{i}  = o(io).V{i}(J);
                    end
                end
                o(io).Idx             = o(io).Idx(J);
                nEl2 = numel(o(io).Idx);
                if nEl1~=nEl2
                    warning(mesId,'%s: %s type = %s, name = %s, numel cells %d > %d (%d vertices removed)',...
                        mfilename,class(o),o(io).type,o(io).name,nEl1,nEl2,nEl1-nEl2);
                end
                if isempty(o(io).Idx)
                    warning(mesId,'%s: %s type = %s, name = %s is empty',...
                              mfilename,class(o),o(io).type,o(io).name);
                    warning(mesId,'off');
                    o(io) = [];
                end
            end
            warning('off',mesId);
        end
        
        function label3(o,varargin)
            %LABEL3: -- place a 3D name label next to the object.
            % Options acctepted by text() can be used. 
            
            for io=numel(o):-1:1
                if ~isempty(varargin)
                    if ~ischar(varargin{1})
                        error('%s: arguments must be [text] followed by prop value pairs compatible with the text function',...
                            mfilename);                    
                    end
                    if isLineSpec(varargin{1})
                        varargin = [o(io).name,varargin{:}];
                    end
                else
                    varargin={o(io).name};
                end
                for i=[1 numel(o(io).P)]
                    text(o(io).P(i).xm,o(io).P(i).ym,o(io).P(i).zm,varargin{:});
                end
            end
        end
        
        function q = Qacross(o,varargin)
            %QACROSS computes the flux [L2/t] acros the object line from left to
            %right looking from its start to its end
            %
            % The total flux across each cell and perpendicular the the
            % line crossing the cell is computed.
            % The results are also integrate over the cells of each line
            % segment the cells take part in.
            % The direction of the flow with more start and end points in a
            % cell, is the vector defined by the first and last point of the
            % lines sections crossing the cell.
            %
            % USAGE:
            %   lineObj = lineObj.qAcross(B(i),gr[,axis,scale],Ilay)
            %   if i is not used, B(end) is used by default
            %   if axis and scale are given, the flows will be plotted as
            %   short lines.
            %
            % TO 131006
            
            [spec     ,varargin] = getWord(varargin,'spec');
            [gr       ,varargin] = getType(varargin,'gridObj',[]);
            [ax       ,varargin] = getType(varargin,'axis',[]);
            [B        ,varargin] = getType(varargin,'struct','');
            [dimension,varargin] = getType(varargin,'char','m3/d');
            
            if isempty(gr)
                error('Needs the grid in the input list');
            end
            
            if ~isempty(ax)
                [scale,varargin] = getType(varargin,'double',1);
                S = dbstack;
                fprintf('%s: will use %g as scale for plot of specific discharge vectors',S(end).name);
            end
            [Ilay ,varargin] = getType(varargin,'double',[]);
            
            if isempty(B)
                error('%s: requires Budget struct and grid',mfilename);
            end

            Nlay = size(B(end).term{1},3);
            if isempty(Ilay)
                Ilay =1:Nlay;
            else
                Ilay(Ilay>Nlay)=Nlay;
                Ilay(Ilay<1   )=1;
                Ilay = unique(Ilay);
            end
            
            for it=numel(B):-1:1
                FRF =  B(it).term{strmatchi('flowR',B(it).label)};
                FFF = -B(it).term{strmatchi('flowF',B(it).label)};
                FRF = sum([FRF(:,1,Ilay), FRF(:,1:end-1,Ilay), FRF(:,end-1,Ilay)],3);
                FFF = sum([FFF(1,:,Ilay); FFF(1:end-1,:,Ilay); FFF(end-1,:,Ilay)],3);

                for io=numel(o):-1:1
                    xp  = [o(io).P.xm];
                    yp  = [o(io).P.ym];
                    dx  = NaN(size(o(io).P));
                    dy  = NaN(size(o(io).P));

                    for ip = numel(o(io).P):-1:1
                        dx(ip) =  diff(o(io).P(ip).x([end 1]));
                        dy(ip) =  diff(o(io).P(ip).y([end 1]));
                    end

                    ds  = sqrt(dx.^2+dy.^2);

                    % Rotate left over 90 degrees
                    dum=   dx;
                    dx =  -dy;
                    dy =   dum;                
                    
                    % Unit vectors, perpendicular to lineObj
                    du  = dx./ds;
                    dv  = dy./ds;

                    q=[];
                    for ip=numel(o(io).P):-1:1
                        
                        for iLay = Nlay:-1:1
                            % specif discharges vectors [L2/T]
                            qx = interp2(gr.xGr,gr.yc,FRF,xp,yp)./gr.DX([o(io).P.idx]);
                            qy = interp2(gr.xc,gr.yGr,FFF,xp,yp)./gr.DY([o(io).P.idx]);
                        end

        %            Use dot product to compute flow perpendicular to line piece in cell
                        q(ip)  = dot([qx(ip); qy(ip)],...
                                     [du(ip); dv(ip)]);
                        if ~isempty(ax)
                            plot(ax,xp(ip)+scale*[0 qx(ip)],...
                                    yp(ip)+scale*[0 qy(ip)],varargin{:});
                        end
                    end
                    o(io).UserData.Qacross(it) = sum(q.*[o(io).P.L]);
                end
            end
            % Handle dimension text depending on object type and request
            % for total flow [L3/T] or object class specific flow, i.e.
            % m3/d for pointObj, m2/d for lineObj and m/d for areaObj
            dims = {
                'm /d',1
                'm /h',1/24
                'M /y',365.24/1e6
                'M /a',365.24/1e6
                'm /y',365.24/1e6
                'm /a',325.24/1e6
                'l /s',1/86.4};
            
            k = find(ismember(dims(:,1),[dimension(1) ' /' dimension(end)]),1,'first');
            if isempty(k)
                error(['Unknown discharge dimension <<%s>>, use oneof\n',...
                            'm3/d m3/h Mm3/a Mm3/y m3/a m3/y l/s'],dimension);
            end

            f = dims{k,2};
            
            if k==size(dims,1)
                if spec
                    dimStr = 'l/s/m';
                else
                    dimStr = 'l/s';
                end
            else
                dimStr = dims{k};
                if spec
                    switch class(o)
                        case 'pointObj', dimStr(2)='3';
                        case 'lineObj' , dimStr(2)='2';
                        case 'area2Obj', dimStr(2)='1';
                    end
                else
                    dimStr(2)='3';
                end
            end
                        
            msgId =  'warning:printQ';
            warning('on',msgId);
            
            t = [B.time];
            
            % Print it
            fprintf('\nQacross sum layers( %s), t [d] -->',sprintf('%d ',Ilay)); 
            fprintf(' %12g',t);
            fprintf('\n');

            for io=1:numel(o)
                
                fldNm = 'Qacross';

                if ~isfield(o(io).UserData,fldNm)
                    warning(msgId,'No fieldName <<%s>> in UserData of %s %s',type_,class(o(io)),o(io).name);
                    continue;
                end
                fprintf('%s %s %4s,',class(o(io)),o(io).name, dimStr);
                if spec
                    fprintf(' %12g',f*sum(o(io).UserData.(fldNm),1)/sum([o(io).P.L]));
                else
                    fprintf(' %12g',f*sum(o(io).UserData.(fldNm),1));
                end
                fprintf('\n\n');
            end
            warning('off',msgId);
        end
    end 
end

classdef area2Obj < lineObj
%AREA2OBJ: generate area2 objects form a spreadsheet.
% area2Obj's are used for area type boundary conditions like lakes pumping
% areas etc.
% area2Obj is a subclass of lineObj which is a subclass of pointObj
%
% USAGE:
%       var = area2Obj([gr,] basename,areaSheetNm,type[options[,plane[,plane[,plane]]]]);
%
% EXAMPLES:
%       riv = area2Obj([gr,] basename,areaSheetNm,gr,'RIV',DEM,[],DEM)
%       drn = area2Obj([gr,] basename,areaSheetNm,gr,'DRN',DEM}
%       ghb = area2Obj([gr,] basename,areaSheetNm,gr,'GHB',DEM}
%       chd = area2Obj([gr,] basename,areaSheetNm,gr,'CHD',DEM,DEM)
%       flux= area2Obj([gr,] basename,areaSheetNm,gr,'WEL',MultiplierArray)
%
% see extended explanation with pointObj (help pointObj)
%    
% SEE ALSO wellObj pintObj lineObj area2Obj
%
% TO 130818
    
    properties
        cellAreaVals  % interpolated values inside polygon
                      % size is   [NcellsInPolygon,size(vertex,o(io).iColx)]
                      % as cellLineVals and vertex.
    end
    properties (Dependent = true)
        perimeter
        area
    end
    methods
        function o=area2Obj(varargin)
        %%AREA2OBj: constructor of area2Obj
        % USAGE:
        %       areas= area2Obj(basename,areaSheetNm,[gr,]type[,values])
        %
        % EXAMPLES:
        %       riv = area2Obj(basename,areaSheetNm,gr,'RIV',DEM,[],DEM)
        %       drn = area2Obj(basename,areaSheetNm,gr,'DRN',DEM}
        %       ghb = area2Obj(basename,areaSheetNm,gr,'GHB',DEM}
        %       chd = area2Obj(basename,areaSheetNm,gr,'CHD',DEM,DEM)
        %       flux= area2Obj(basename,areaSheetNm,gr,'WEL',MULTARRAY)
        %
        % TO 130813
            
            o = o@lineObj(varargin{:});
            
            if nargin<2, return; end
            
            %% Process input in the same manner as in super class to make sure
            %  we can reprocess some of the argument in a subclass-specific way
            [~ ,varargin] = getProp(varargin,'active','');
            [~ ,varargin] = getProp(varargin,'conc',[]);
            [~ ,varargin] = getProp(varargin,'struct',{});

            [gr,varargin] = getType(varargin,'gridObj',[]);
            [~ ,varargin] = getNext(varargin,'char',[]);
            [~ ,varargin] = getNext(varargin,'char',[]);
            [~ ,varargin] = getNext(varargin,'char',[]);

            [colHdr_  ,varargin] = getNext(varargin,'char','');
            if ~isempty(colHdr_)
                [~ ,varargin] = getNext(varargin,{'char','double','cell'},[]);
            end
            
            if isempty(gr),return; end
            
            %% ======= Put areas into grid ========

            for io = numel(o):-1:1           
                
                if isempty(o); break; end
                
                %% Interpolate contour data to internal area
                   [o(io).cellAreaVals,o(io).Idx] = ...
                       gr.interpLaplacian([o(io).iColx,o(io).iColy],o(io).cellLineVals,[o(io).P.idx]);
                
                % add missing columns (id,x,y) using xm and ym
                o(io).A = gr.AREA3(o(io).Idx(:,1));
                
                %% processing further command-line arguments.
                % These are assumed to be planes of size [Ny,Nx]
                % or a scalar.
                
                for iPlane = 1:min(o(io).narg,numel(varargin))
                    if isempty(varargin{iPlane})
                        continue;
                    end
                    % Maybe use a scalar instead of a plane. TO 141103
                    if isscalar(varargin{iPlane})
                        o(io).V{iPlane} = varargin{iPlane};
                        continue;
                    end
                    if ~all(size(varargin{iPlane}(:,:,1))==[gr.Ny,gr.Nx])
                        error('%s: argument must be of size Ny,Nx = [%d,%d]',...
                            mfilename,gr.Ny,gr.Nx);
                    end
                    
                    o(io).V{iPlane} = varargin{iPlane}(o(io).Idx);
%                     switch o(io).type
%                         case 'WEL'
%                             if iPlane == 1
%                                 o(io).V{iPlane} = o(io).V{iPlane}.*o(io).A;
%                             end
%                         case {'DRN','GHB','RIV'}
%                             if iPlane == 2
%                                 o(io).V{iPlane} = o(io).V{iPlane}.*o(io).A;
%                             end
%                         case 'CHD'
%                             % skip
%                         otherwise
%                     end
                end
                
            end
        end
               
        function perimeter = get.perimeter(o)
            %PERIMETER-- computes the perimeter of the area object
            perimeter = ...
                sum(sqrt(diff(o.vertex(:,o.iColx)).^2 + diff(o.vertex(:,o.iColy)).^2));                    
        end
        
        function area = get.area(o) 
            %AREA: computes surface area of a polyline considered as polygon
            % computation based on outer products of vectors spanned
            % between point in polygond and each line segment.
            xv=o.vertex(2:end-1,o.iColx);
            yv=o.vertex(2:end-1,o.iColy);

            xp = o.vertex(1,o.iColx);
            yp = o.vertex(1,o.iColy);

            dx = xv-xp;
            dy = yv-yp;

            dA = NaN(3,numel(dx)-1);

            for i = size(dA,2):-1:1
                dA(:,i) =  cross([dx(i); dy(i); 0], [dx(i+1); dy(i+1); 0]);
            end

            area = abs(sum(dA(3,:))/2);
        end
 
        function o = plotPnts(o,varargin)
            %PLOTPNTS: plots markers at center of cells in areaObj
            % TO 130816
            [ax,varargin] = getNext(varargin,'axis',gca);
            plot(ax,gr.XM(o(io).Idx),gr.YM(o(io).Idx),'o',varargin{:})
        end
        
        function o = check(o,IBOUND)
            %CHECK: remove cells that conincide with IBOUND==0
            % TO 130816
            mesId = [class(o) ':pointsAtZeroIBOUND'];
            warning('on',mesId);
            for io=numel(o):-1:1
                if nargin<2 || ~all(size(IBOUND)==o(io).grSize(1:numel(size(IBOUND))))
                    error('%s: need IBOUND in call and its size must be %s',...
                        mfilename,sprintf(' %d',o(iL).size));
                end
                
                % first check the cells on the boundary of the area
                nEl1 = numel(o(io).P);
                J = find(IBOUND([o(io).P.idx])~=0);
                o(io).cellLineVals    = o(io).cellLineVals(J,:);
                o(io).P               = o(io).P(J);
                o(io).sCell           = o(io).sCell(J);
                nEl2 = numel(o(io).Idx);
                if nEl1~=nEl2
                    warning(mesId,'%s: %s type = %s, name = %s, numel cells %d > %d (%d vertices removed)',...
                        mfilename,class(o),o(io).type,o(io).name,nEl1,nEl2,nEl1-nEl2);
                end
                if isempty(o(io).P)
                    warning(mesId,'%s: %s type = %s, name = %s is empty',...
                              mfilename,class(o),o(io).type,o(io).name);
                    warning(mesId,'off');
                    o(io) = [];
                    continue;
                end
                
                % then check the cells within the area
                nEl1 = numel(o(io).Idx);
                J = find(IBOUND(o(io).Idx)~=0);
                o(io).cellAreaVals    = o(io).cellAreaVals(J,:);
                for i=1:numel(o(io).V)
                    if numel(o(io).V{i})>1
                        o(io).V{i} = o(io).V{i}(J);
                    end
                end
                o(io).A = o(io).A(J);
                o(io).Idx             = o(io).Idx(J);
                nEl2 = numel(o(io).Idx);
                if nEl1~=nEl2
                    warning(mesId,'%s: %s type = %s, name = %s, numel cells %d > %d (%d removed)',...
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
        
        function  o = fill(o,varargin)
            %FILL fills area2Obj with specified or default color
            %
            % EXAMPLE:
            %         area2Obj = area2Obj.fill();
            %         area2Obj = area2Obj.fill('w');
            %         area2Obj = area2Obj.fill(ax,'r','faceAlpha',0.3);
            %
            % TO 130825
            
            [ax,varargin] = getNext(varargin,'axis',gca);
            
            % is lineSpec specified ?? If so use it for all objects
            colorSpecSpecified = ~isempty(varargin) && isColor(varargin{1});
            if colorSpecSpecified
                color = varargin{1};
                varargin(1) = [];
            end
            
            for io=1:numel(o)
                if ~colorSpecSpecified
                    color =o(io).EdgeColor;
                end
                patch(o(io).vertex(:,o(io).iColx),o(io).vertex(:,o(io).iColy),color,varargin{:},'parent',ax);
            end
        end

    end
end

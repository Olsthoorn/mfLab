classdef areaObj < esriShapeObj
    properties
        x,y
        nr, id, izone
        name, shortname, longname
        code
        iz     % layers in which area object is active
        idx,   % indices of areaObj in first layer
        LRC    % Lay Row Col of this waterbody in mesh (redundant)
        UserData
    end
    properties(Dependent = true)
        box
        numPoints
        numParts
        zRange
        mRange
        xc
        yc
        perimeter
        area
        edge
    end
    methods
        function o=areaObj(varargin)

            if nargin<1; return; end
            
            % Has the grid been passed to this constructor?
            % if so, get it and make sure toGrid is run
            for i=1:numel(varargin)
                if strcmpi(class(varargin{i}),'gridObj')
                    gr = varargin{i};
                    varargin(i)=[];
                end
            end
            
            switch numel(varargin)
                case 3  % USAGE areaObj(basename,sheetNmProps,sheetNmCoords) 
                    if ~ischar(varargin{1})
                        error('%s: first argument must be a file name not of class <<%s>>',class(varargin{1}));
                    end
                    
                    
                    [STATUS,SHEETS] = xlsfinfo(varargin{1});
                    if ~isempty(STATUS)
                        if ~strmatchi(varargin{2},SHEETS)
                            error('%s: Can''t fine sheet <<%s>> in workbook <<%s>>',...
                                mfilename,varargin{2},varargin{1});
                        end
                        if ~strmatchi(varargin{3},SHEETS)
                            error('%s: Can''t find sheet <<%s>> in workbook <<%s>>',...
                                mfilename,varargin{3},varargin{1});
                        end
                        
                        xlsName  = varargin{1};
                        sheetNm1 = varargin{2};
                        sheetNm2 = varargin{3};
                        %% get data from spreadsheet
                        
                        % fieldnames AreaObj
                        members = fieldnames(o);
                        
                        % get info from sheetNm1 into table
                        table = getTable(xlsName,sheetNm1,members,'horizontal');
                        
                        % found fields
                        field = fieldnames(table);

                        %% get coordinates from second sheet name
                        [PNm,Coords] = getExcelData(xlsName,sheetNm2,'Vertical');

                        % get rows starting with label x (ie. x1 x2 x3 etc)
                        Ix = strmatchi('x',PNm);
                        
                        %% assert there are coordinates for each areaObj
                        if numel(Ix) ~= numel(table)
                            error(['%s: The number of lines with x-coordinates <<%d>> in workbook <<%s> worksheet <<%s>>',...
                                'does not match the number of areadObj <<%s>> in workbook <<%s>> worksheet <<%s>>'],...
                                mfilename,numel(Ix),xlsName,sheetNm1,numel(table),xlsName,sheetNm2);
                        end
                        
                        % for all areaObjects
                        for io = numel(table):-1:1
                            
                            % fill the fields
                            for j=1:numel(field)
                                o(io).(field{j}) = table(io).(field{j});
                            end
                        
                            % fill the coordiinates into the areaObj
                            j = min(find(~isnan(Coords(Ix(io)  ,:)),1,'last'),...
                                    find(~isnan(Coords(Ix(io)+1,:)),1,'last'));
                            if ~isempty(j),
                                o(io).x = Coords(Ix(io)  ,1:j);
                                o(io).y = Coords(Ix(io)+1,1:j);
                            else
                                o(io).x = Coords(Ix(io)  ,:);
                                o(io).y = Coords(Ix(io)+1,:);
                            end
                            
                            % set some fields
                            o(io).type      = 'area';
                            o(io).typeName  = 'area';
                            o(io).parts     = [];
                            o(io).points    = [o(io).x(:) o(io).y(:)];
                            o(io).zArray    = [];
                            o(io).mArray    = [];
                            o(io).created   = now;
                            o(io).filename  = xlsName;
                        end
                    end
                otherwise
            end
            
            for io = numel(o):-1:1
              o(io).izone = o(io).nr; o(io).id = o(io).nr;
            end
            
            if exist('gr','var')
                o = o.toGrid(gr);
            end
        end
        
        function numPoints = get.numPoints(o), numPoints = size(o.points,1); end
        function numParts  = get.numParts(o) , numParts  = size(o.parts,1); end
        function zRange    = get.zRange(o),    zRange    = [min(o.zArray(:)) max(o.zArray(:))]; end
        function mRange    = get.mRange(o),    mRange    = [min(o.mArray(:)) max(o.mArray(:))]; end
        function xc  = get.xc(o), xc=mean(o.points(:,1)); end
        function yc  = get.yc(o), yc=mean(o.points(:,2)); end
        function box = get.box(o)
            box=[min(o.points(:,1)),min(o.points(:,2)),max(o.points(:,1)),max(o.points(:,2))];
        end
        function edge= get.edge(o)
            N = o.numPoints;
            edge(N).x = [o.points(end,1) o.points(1,1)];
            edge(N).y = [o.points(end,2) o.points(1,2)];
            for i=(N-1):-1:1
                edge(i).x = [o.points(i,1) o.points(i+1,1)];
                edge(i).y = [o.points(i,2) o.points(i+1,2)];
                edge(i).L = sqrt(diff(edge(i).x).^2 + diff(edge(i).y).^2);
            end
        end
        function perimeter = get.perimeter(o)
            % peremiter = areaObj.peremiter(); -- peremiter of area object
            perimeter= sum([o.edge.L]);
        end
        function area=get.area(o)
            % area = areaObj.area(); -- exact area of area object
            dx=o.points(:,1)-o.xc;
            dy=o.points(:,2)-o.yc;
            area=sum(cross([dx,            dy,            zeros(size(dx))],...
                           [dx([2:end 1]), dy([2:end 1]), zeros(size(dx))]))/2;
            area=abs(area(end));
        end

        function o=rotate(o,xc,yc,alfa)
            % rotate object around point xc,yc over alfa counter clock wise
            [o.points(:,1) o.points(:,2)]=mf_rotate(o.points(:,1)-xc,o.points(:,2)-yc,0,0,alfa);
        end
        function plot(o,varargin)
            % areaObj.plot(clrs) -- plot areaObjects
            [ax  ,varargin] = getProp(varargin,'axis',gca);
            [ax  ,varargin] = getNext(varargin,'axis',ax);
            [clrs,varargin] = getNext(varargin,'color','brgkmc');
            
            for io = 1:numel(o)
                patch(o(io).points(:,1),o(io).points(:,2),mf_color(io,clrs),'parent',ax,varargin{:});
            end
        end
        
        function whoami(o); fprintf('I am a <<%s>>\n',class(o)); end

        function o= toGrid(o,gr)
            % o= toGrid(o,gr,izone); -- put areaObj into the grid (layer 1)
            warning('on','matLab:areaObj:areaOutsideGrid');
            for i = numel(o):-1:1
                o(i).iz    = 1;
                o(i).idx   = find(inpolygon(gr.XM(:,:,1),gr.YM(:,:,1),o(i).points(:,1),o(i).points(:,2)));
                if isempty(o(i).idx)
                    warning('matLab:areaObj:areaOutsideGrid',...
                        ['%s: area <<%s>> nr <<%d>> izone <<%d>> lies entirely outside the grid.\n',...
                        'Could it be that the area is so small that there is no cell center falling in it?\n'],...
                        mfilename,o(i).name,o(i).nr,o(i).izone);
                end
                o(i).idx = o(i).idx(:);
                o(i).LRC = cellIndices(o(i).idx,gr.size,'LRC');
            end
        end
        
        function o = setLRC(o,gr), o.LRC=cellIndices(o.Idx,gr.size,'LRC'); end
        
        function ZONE = ZONE(o,gr)
            % ZONE = areaObj.ZONE(gr); -- generates a zone array where
            % using the izone of the areaObj as zone numbers (defalt areaObj.izone == areaObj.nr)
            ZONE = zeros(gr.Ny,gr.Nx);
            for io=1:numel(o)
                ZONE(o(io).idx) = o(io).nr;
            end
        end
        
    end
end
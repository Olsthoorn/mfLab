classdef waterBodyObj < areaObj
    properties
        stage   % water level
        Idx     % indices of areaObj down to layer iz(end)
        c       % bottom entry resistance
        z_bot   % bottom of the area object (of the lake for instance, const)
        gauge   % gauge       
    end
    properties (Dependent = true)
        vcont % vertical conductance = d/k = 1/c
    end
    
    methods
        function o=waterBodyObj(varargin)
                        
            o = o@areaObj(varargin{:});
            
            if nargin==0; return; end
            
            % Has the grid been passed to this constructor?
            % if so, get it and make sure toGrid is run
            for i=numel(varargin):-1:1
                if strcmpi('gridObj',class(varargin{i}))
                    gr = varargin{i};
                    hasGrid = true;
                    break;
                else
                    hasGrid = false;
                end
            end

            for i=1:numel(o)
                if isfield(o(i).UserData,'h'     ),o(i).stage   = o(i).UserData.h; end
                if isfield(o(i).UserData,'z_bot' ),o(i).z_bot   = o(i).UserData.z_bot; end
                if isfield(o(i).UserData,'xGauge'),o(i).gauge.x = o(i).UserData.xGauge; end
                if isfield(o(i).UserData,'yGauge'),o(i).gauge.y = o(i).UserData.yGauge; end

                o(i).izone = o(i).nr;

                if hasGrid, o(i) = o(i).toGrid(gr); end
            end
        end
    
        function vcont = get.vcont(o)
            % vertical conductivity of layer i.e. d/kv [1/T] or 1/c (resistance [T])
            vcont = 1./o.c;
        end

        function o = setKML(o,kmlfile)            
            xy = wgs2rd(kmlpath(kmlfile));
            o.x = xy(:,1);
            o.y = xy(:,2);
        end

        function o = setXY(o,varargin)
            % waterbodyObj.setXY([x y]); set coordinates of circumference
            % of waterbody
            switch numel(varargin)
                case 1
                    o.x = varargin{1}(:,1);
                    o.y = varargin{1}(:,2);
                otherwise
                    o.x = varargin{1};
                    o.y = varargin{2};
            end
        end
        
        function o = setDEM(o,gr,varargin)
            % waterbody = waterBody.setDem(gr,aDem);
            if gr.Nx ~= size(varargin{1},2)
                error('%s: size of DEM along columns does not mach size of model',mfilename);
            end
            if gr.Ny ~= size(varargin{1},1)
                error('%s: size of DEM along rows does not match size of model',mfilename);
            end
            o.DEM = varargin{1};
            
            IN = inpolygon(gr.XM,gr,YM,o.x,o.y);
            o.DEM = o.DEM & IN;
        end

        function zarray = zoneArray(o,gr)
            % ZONE = areaObj.ZONE(gr); -- generates a zone array where
            % using the izone of the areaObj as zone numbers (defalt areaObj.izone == areaObj.nr)
            zarray = gr.const(0);
            for io=1:numel(o)
                zarray(o(io).Idx) = o(io).nr;
            end
        end

        function o= toGrid(o,gr)
            % o= toGrid(o,gr,izone); -- put areaObj into the grid
            warning('on','matLab:areaObj:areaOutsideGrid');
            for i = numel(o):-1:1
                o(i).iz    = xyzindex(o(i).z_bot,gr.zGr);
                o(i).idx   = find(inpolygon(gr.XM(:,:,1),gr.YM(:,:,1),o(i).x,o(i).y));
                if isempty(o(i).idx)
                    warning('matLab:areaObj:areaOutsideGrid',...
                        ['%s: area <<%s>> nr <<%d>> izone <<%d>> lies entirely outside the grid.\n',...
                        'Could it be that the area is so small that there is no cell center falling in it?\n'],...
                        mfilename,o(i).name,o(i).nr,o(i).izone);
                end
                o(i).idx   = o(i).idx(:);
                n       = numel(o(i).idx);

                if ~isempty(o(i).idx)
                    o(i).Idx(n*o(i).iz,1)=NaN;
                    for iz_=1:o(i).iz
                        o(i).Idx((iz_-1)*n+(1:n),1) = gr.Ny*gr.Nx*(iz_-1)+o(i).idx(:);
                    end
                    fprintf('waterbody %s nr(%d) in mesh\n',o(i).shortname,o(i).izone);
                end
                o(i).LRC = cellIndices(o(i).Idx,gr.size,'LRC');
            end
        end
        
        function h=head(o,Piez,time)
            k=strmatchi(o.gauge,{Piez.name},'exact');
            if k(1)>0
                h=interp1([-Inf; Piez(k).t(:); Inf] , Piez(k).h([1 1:end end] ),time);
            else
                error('Can''t find Piezometer <<%s>> for waterbody <<%s>>,Piez(%d)',Piez.name,o.name,o.gauge);
            end
        end

    end
end
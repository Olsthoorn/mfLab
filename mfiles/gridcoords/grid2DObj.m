classdef grid2DObj
    %GRID2DOBJ -- 2D gridObj for use with finite difference models
    % USAGE: gr = grid2DObj(xGr,yGr[,Z],'AXIAL'[,true|{false}])
    %  xGr = grid line coordinates of columns
    %  yGr = grid line coordinates of rows
    %  Z   = vector or 3D array of tops and bottoms of layer
    %        may be vector, [1 0] is assumed when omitted.
    % 'AXIAL',true|false indicates whether the grid is for an
    %        axially symmetric model or not, default = false
    % 'AXIAL' without true or false, then true is assumed
    properties
        % stored properties
        xGr,yGr,zGr,Z
        AXIAL = false;
        LAYCON  = false;  % LAYCON in z-direction, horizontal model
        LAYCON2 = false;  % LAYCON in y-direction, vertical   model 
    end
    properties (Dependent=true)
        % properties computed when used
                        rGr
        Nx,Ny,Nz,Nod,   Nr
        XGr,YGr,        RGr
        dx,dy,dz,       dr
        dX,dY,dZ,       dR
        xm,ym,zm,       rm
        Xm,Ym,Zm,       Rm
        xp,yp,          rp
        Xp,Yp,          Rp
        Area,AREA
        Vol,VOL
    end
    methods
        function o = grid2DObj(varargin)
            %GRID2DOBJ -- constructor
            % USAGE: gr = grid2DObj(xGr,yGr[,Z],'AXIAL'[,true|{false}])
            %  xGr = grid line coordinates of columns
            %  yGr = grid line coordinates of rows
            %  Z   = vector or 3D array of tops and bottoms of layer
            %        may be vector, [1 0] is assumed when omitted.
            % 'AXIAL',true|false indicates whether the grid is for an
            %        axially symmetric model or not, default = false
            % 'AXIAL' without true or false, assumes true
            % 'LAYCON' true|false  aquifer is watertable if LAYCON == true
            % 'LAYCON' without true or false, assumes true
            % 'LAYCON2' same but for cross section models (where y is in
            %           fact the vertical z direction.
            % TO 140401
            
            if nargin==0, return; end
            
            [o.AXIAL  ,varargin]  = getProp(varargin,'AXIAL',false);
            if ~o.AXIAL, [o.AXIAL ,varargin]=getWord(varargin,'AXIAL'); end

            [o.LAYCON2,varargin]  = getProp(varargin,'LAYCON2',false);            
            if ~o.LAYCON2,[o.LAYCON2,varargin] = getWord(varargin,'LAYCON2'); end
            
            if ~o.LAYCON2
                [o.LAYCON ,varargin]  = getProp(varargin,'LAYCON',false);                        
                if ~o.LAYCON ,[o.LAYCON ,varargin] = getWord(varargin,'LAYCON'); end
            end
            
            if o.AXIAL,   fprintf('Notice: grid has AXIAL = true\n'); end
            if o.LAYCON,  fprintf('Notice: grid has LAYCON = true\n'); end
            if o.LAYCON2, fprintf('Notice: grid has LAYCON2 = true\n'); end
            
            [xGr_,varargin] = getNext(varargin,'double',[]);
            [yGr_,varargin] = getNext(varargin,'double',[]);
            [Z_  ,    ~   ] = getNext(varargin,'double',[1 0]);
                        
            if ~isempty(varargin)
                fprintf('Notice: not all inputs were used !\n');
            end
            
            o.xGr = unique(xGr_(:)');
            o.yGr = unique(yGr_(:) ); o.yGr = o.yGr(end:-1:1);
            
            if nargin<3 || o.AXIAL
                o.Z   = cat(3,o.const(1),o.const(0));
                o.zGr = mean(mean(o.Z,1),2);
            elseif size(Z_,3)>1
                o.Z   = Z_(:,:,1:2);
                o.zGr = mean(mean(o.Z,1),2);
            else
                if isvector(Z_)
                    o.Z = XS(unique(Z_(:)));
                    o.Z = o.Z(end:-1:1);
                end
                o.zGr = mean(mean(o.Z,1),2);
                if size(o.Z,1)~=o.Ny
                    o.Z = bsxfun(@times,ones(o.Ny,1),o.Z);
                end
                if size(o.Z,2)~=o.Nx
                    o.Z = bsxfun(@times,ones(1,o.Nx),o.Z);
                end
            end
            % verify
            if o.Ny<1, error('Need at least 2 distinct yGr coordinates.'); end
            if o.Nx<1, error('Need at least 2 distinct xGr coordinates.'); end
            if o.Nz<1, error('Need at least 2 distinct zGr coordinates.'); end
            if o.Nz>1, error('Use at most 2 z-elevations in 2D models.');  end
        end
        function sz = size(o,dim)
            %SIZE -- size of grid
            % USAGE: sz = gr.size([dim])
            %   dim is optional or 1 or 2
            if nargin<2,  sz = [o.Ny o.Nx]; return; end
            
            if  dim<2, sz = o.Ny;  else  sz = o.Nx;  end
        end
        function A = const(o,a)
            %CONST -- generate constant array of size Ny,Nx filled with a
            % USAGE A = const(a)
            % a is a scalar or first entry of array or vector
            if isscalar(a)
                A = a + zeros(o.Ny,o.Nx);
            elseif size(a,2) == 1
                A = a * ones(1,o.Nx);
            elseif size(a,1) == 1
                A = ones(o.Ny,1) * a;
            else
                error('input must be scalar or a vector');
            end
        end
        function plot(o,varargin)
            %PLOT -- plots the grid/mesh using linespec
            % USAGE: gr.plot([lineSpec[,property,value[,property,value[,...]]]])
            %   lineSpec is legal Matlab line specification for plot (see
            %   plot function in Matlab, e.g 'r'  'b--'  'ko.-')
            %   property value pairs should be legal propery values pairs
            %   for plot function in Matlab (e.g, 'lineWidth',2, ...)
            set(gca,'nextPlot','add');
            for ix=1:numel(o.xGr)
                plot(o.xGr([ix ix]),o.yGr([1 end]),varargin{:});
            end
            for iy=1:numel(o.yGr)
                plot(o.xGr([1 end]),o.yGr([iy iy]),varargin{:});
            end
        end
        function [xr,Ix] = xr(o,x)
            %XR -- relative coordinates left oriented
            % USAGE: [xr,Ix] = gr.xr(x);
            Ix = interp1(o.xGr,1:o.Nx+1,x);
            xr = Ix - floor(Ix);
            Ix =      floor(Ix);
            xr(Ix==o.Nx+1) = 1;
            Ix(Ix==o.Nx+1) = o.Nx;             
        end
        function [xrm,Ix] = xrm(o,x)
            %XRM -- relative coordinages, central oriented
            % USAGE: [xrm,ix] = gr.xrm(x);
            [xrm,Ix] = o.xr(x);
            xrm = xrm-0.5;
        end
        function [yr,Iy] = yr(o,y)
            %YR -- relative y coordinates, top oriented
            % USAGE: [yr,iy] = gr.yr(y);
            Iy = interp1(o.yGr,1:o.Ny+1,y);
            yr = Iy - floor(Iy);
            Iy =      floor(Iy);
            yr(Iy==o.Ny+1) = 1;
            Iy(Iy==o.Ny+1) = o.Ny;
        end
        function [yrm,Iy] = yrm(o,y)
            %YRM -- relative coordinates, center oriented
            % USAGE: [yrm,iy] = gr.yrm(y)
            [yrm,Iy] = o.yr(y);
            yrm = yrm - 0.5;
        end
        function [Idx,Ix,Iy,xr,yr,xrm,yrm] = Idx(o,x,y)
            %IDX -- global index of points x,y
            % USAGE: [Idx,Ix,Iy,xr,yr,xrm,yrm] = Idx(x,y)
            [xr,Ix] = o.xr(x);
            [yr,Iy] = o.yr(y);
            Idx = o.Ny*(Ix-1)+Iy;
            xrm = xr - 0.5;
            yrm = yr - 0.5;
        end
        function Ix = Ix(o,x), Ix = max(1,ceil(interp1(o.xGr,0:o.Nx,x))); end
        function Iy = Iy(o,y), Iy = max(1,ceil(interp1(o.yGr,0:o.Ny,y))); end
        function ux = ux(o,x), [~,ux] = o.xrm(x); end
        function uy = uy(o,y), [~,uy] = o.yrm(y); end
        
        % Dependent variables
        function Ny = get.Ny(o), Ny = numel(o.yGr)-1; end
        function Nx = get.Nx(o), Nx = numel(o.xGr)-1; end
        function Nz = get.Nz(o), Nz = numel(o.zGr)-1; end
        function Nod= get.Nod(o), Nod = o.Ny*o.Nx*o.Nz; end
        function dx = get.dx(o), dx = abs(diff(o.xGr,1,2)); end
        function dy = get.dy(o), dy = abs(diff(o.yGr,1,1)); end
        function dz = get.dz(o), dz = abs(diff(o.zGr,1,3)); end
        function xm = get.xm(o), xm = 0.5*(o.xGr(1:end-1) + o.xGr(2:end)); end
        function ym = get.ym(o), ym = 0.5*(o.yGr(1:end-1) + o.yGr(2:end)); end
        function zm = get.zm(o), zm = 0.5*(o.zGr(1:end-1) + o.zGr(2:end)); end
        function xp = get.xp(o), xp = o.xGr(2:end-1); end
        function yp = get.yp(o), yp = o.yGr;          end
        function XGr = get.XGr(o), XGr = ones(size(o.yGr)) * o.xGr; end
        function YGr = get.YGr(o), YGr = o.yGr * ones(size(o.xGr)); end
        function dX  = get.dX(o), dX = ones(o.Ny,1) * o.dx; end
        function dY  = get.dY(o), dY = o.dy * ones(1,o.Nx); end
        function dZ  = get.dZ(o), dZ = abs(diff(o.Z  ,1,3)); end
        function Xm  = get.Xm(o), Xm = ones(o.Ny,1) * o.xm; end
        function Ym  = get.Ym(o), Ym = o.ym * ones(1,o.Nx); end
        function Zm  = get.Zm(o), Zm = 0.5*(o.Z(:,:,1:end-1)+o.Z(:,:,2:end)); end
        function Xp  = get.Xp(o), Xp = ones(size(o.yp)) * o.xp; end
        function Yp  = get.Yp(o), Yp = o.yp * ones(size(o.xp)); end
        function Area = get.Area(o)
            if ~o.AXIAL, Area = o.dy * o.dx;
            else         Area = pi.*(o.xGr(2:end).^2-o.xGr(1:end-1).^2);
            end
        end
        function AREA = get.AREA(o), AREA = sum(o.Area(:)); end
        function Vol = get.Vol(o)
            if ~o.AXIAL, Vol = bsxfun(@times,o.Area,o.dZ);
            else         Vol = bsxfun(@times,o.Area,o.dy);
            end
        end
        function VOL  = get.VOL(o),  VOL = sum(o.Vol(:));  end

        %aliases for x-variables in axially symmetric cases
        function Nr  = get.Nr(o),  Nr  = o.Nx;  end
        function rGr = get.rGr(o), rGr = o.xGr; end
        function RGr = get.RGr(o), RGr = o.XGr; end
        function rm  = get.rm(o),  rm  = o.xm;  end
        function Rm  = get.Rm(o),  Rm  = o.Xm;  end
        function dr  = get.dr(o),  dr  = o.dx;  end
        function dR  = get.dR(o),  dR  = o.dX;  end
        function rp  = get.rp(o),  rp  = o.xp;  end
        function Rp  = get.Rp(o),  Rp  = o.Xp;  end
    end
end
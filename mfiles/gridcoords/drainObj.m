classdef drainObj < hydObj
%DRAIN class def for hydObj
%
%  ToDo: to be worked out. Now preliminarily used by Drain16
%
% TO 123410
    
    properties
        xheel; yheel; zheel % Drain outleet position
        z_bot   % bottom of the drain
        c_entry % drain entry resistance
        t; h    % drain head time series        
        rw      % drain radius
        Rw
        c0; cL; lambda; dphi %drain resistance parameters
        L; dL; I;
        %LRC;
        %ix;
        %iy;
        %iz;
        xm;
        ym;
        zm  % Drain position
        kDrLong % cond in longitudinal direction (along drain)
        HKx     % drain cell cond along            drain horizontal    (x direction)
        HKy     % drain cell cond perpendicular to drain horizontal
        VK      % drain cell cond perpendicular to drain vertical
        STRTHD  % drain start head
        drainhd % dain iniital head
        c       % drain resistance
    end
    
    methods
        function o = drainObj(shape,varargin)
            %DRAINOBJ drain object, to be worked out.
            %
            % ToDo: to be worked out. Now preliminarily used by Drain16

            if nargin==0, return; end
            
            if nargin>1
                o.name=shape;
                ni=numel(varargin);
                if ni>=2
                    o.name   = shape;
                    o.shortnm= shape;
                    o.x      = varargin{1};
                    o.y      = varargin{2};
                    o.z      = varargin{3};
                end
                if ni>=4
                    o.code   = varargin{4};
                end
                if ni>=5
                    o.xheel  = varargin{5};
                    o.yheel  = varargin{6};
                    o.zheel  = varargin{7};
                end
                return;
            end

            switch nargin
                case 0, return;
                case 1, % interprete shape structure array
                    o(numel(shape),1)=drainObj();
                    for i=numel(shape):-1:1
                        o(i).type   = 'drain';
                        o(i).name   = shape(i).NAAM;
                        o(i).x      = shape(i).x{:};
                        o(i).y      = shape(i).y{:};
                        o(i).code   = shape(i).CODE;
                    end
                otherwise
                    error('%s constructor: Wrong number of arguments.',class(o));
            end
        end
        function o = rotate(o,xc,yc,alfa)
            [o.x ,o.y ]=mf_rotate(o.x -xc,o.y -yc,0,0,alfa);
            [o.xheel,o.yheel]=mf_rotate(o.xheel-xc,o.yheel-yc,0,0,alfa);
        end
        function h = head(o,Piez,time)
            h=...
                interp1([-Inf; Piez.t(:); Inf],Piez.h([1 (1:end) end]),time) *ones(size(o.I));
                  %  + o.dphi*(exp(o.L/max(o.L)).^o.n-1)/(exp(1)-1);   Do this later     
        end
        function o = merge(o,grid)
            % Draw line through grid
            PDr=linegrid([o.x(:) o.y(:) o.z(:)],grid.xGr,grid.yGr,grid.zGr);
            o.LRC = [[PDr.iz]' [PDr.iy]' [PDr.ix]'];
            o.I  =   [PDr.I];
            o.ix =   [PDr.ix];
            o.iy =   [PDr.iy];
            o.iz =   [PDr.iz];
            o.xm =   [PDr.xm];
            o.ym =   [PDr.ym];
            o.zm =   [PDr.zm];
            o.L=sqrt((o.xm-o.xheel).^2 ...
                      +(o.ym-o.yheel).^2 ...
                      +(o.zm-o.zheel).^2);
            
            % Computing the conductivity of drain cells so that the specified entry resistance
            % is met exactly, see theory

            % the conductivity that incorporates the drain ressitance has been derived
            % analytically and reads k = pi R/c / (DY/DZ+DZ/DY)
            
            % But the resistance is a function of distance from heel fixed
            % by c0(heel), cL(toe) and lambda [m] to scale the curve's
            % curvature
            
            c1=(o.c0-o.cL*exp(-max(o.L)/o.lambda))/(1-exp(-max(o.L)/o.lambda));
            o.c=c1+(o.cL-c1)*exp(-(max(o.L)-o.L)/o.lambda);
            
            o.VK =(pi*o.rw./o.c)./(grid.DY(o.I)./grid.DZ(o.I)+grid.DZ(o.I)./grid.DY(o.I));
            o.HKy=o.VK;
            o.HKx=o.kDrLong;
        end
    end
    
end


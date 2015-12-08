classdef hydobj
    %hydobj general hydrological object
    %   Detailed explanation goes here
    properties
        type
        name
        shortname
        longname
        gauge
        code
        nr
        x
        y
        tstart; tend
    end
    properties (Dependent=true)
        box
        XC
        YC
        perimeter
        area
    end
    methods
        function o=hydobj()
            switch nargin
                case 0, return;
                otherwise
                    error('%s constructor, wrong number of arguments(%d)',nargin);
            end        
        end
        
        function o=setXY(o,x,y)
            if nargin==2
                if ~size(x,2)==2,
                        error('%s constructor, dimension of x must be 2 column array in this case',class(o));
                else
                    y=x(:,2); x=x(:,1);
                end
            end
            o.x=x(:)';
            o.y=y(:)';
        end
        
        function o=shift(o,dx,dy)
            o.x=o.x+dx;
            o.y=o.y+dy;
        end
        
        function XC=get.XC(o), XC=mean(o.x); end
        function YC=get.YC(o), YC=mean(o.y); end
        
        function box=get.box(o)
           box=[min(o.x),min(o.y),max(o.x),max(o.y)];
        end
        
        function perimeter=get.perimeter(o)
                perimeter=sqrt(sum(diff(o.x).^2+diff(o.y).^2));
        end
        function area=get.area(o)
            dx=o.x(:)-o.XC; dy=o.y(:)-o.YC;
            area=sum(cross([dx             dy           zeros(size(dx))],...
                           [dx([2:end 1]) dy([2:end 1]) zeros(size(dx))]))/2;
            area=abs(area(end));
        end
        
        function o=rotate(o,xc,yc,alfa)
            [o.x ,o.y ]=rotate(o.x-xc,o.y-yc,0,0,alfa);
        end
        
        function plot(o,cl,lw,type)
            if nargin<4, type=o.type; end
            switch nargin
                case 1, cl='b-'; lw=1;
                case 2, lw=1;
                case 3,
                case 4
                otherwise
                    error('%s.plot, wrong number of arguments',class(o));
            end
            if strcmpi(type,o.type)
                plot(o.x,o.y,cl,'linewidth',lw);
            end
        end
        
        function whoami(o)
            fprintf('I am a <<%s>>\n',class(o));
        end
    end 
end

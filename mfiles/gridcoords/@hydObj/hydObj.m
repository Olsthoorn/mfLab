classdef hydObj
    %HYDOBJ class def for general hydrological object
    %   Detailed explanation goes here
    properties
        nr
        id
        type = 'genearl hydObj'
        name; shortnm; longname;
        code
        x
        y
        z
        ix, iy, iLay, iz, idx
        LRC    % Lay Row Col of this waterbody in mesh (redundant)

        UserData
    end
    properties (Dependent=true)
        box
        xc
        yc
        perimeter
        area
        edge
    end
    methods
        function o=hydObj(nr,x,y,z)
            %HYDOBJ constructor for hydrological object
            %
            switch nargin
                case 0,
                    return
                case 1 % <[x y]>
                    nr  = nr(:);
                    o.x = nr(:,1)';
                    o.y = nr(:,2)';
                case 2 %  <x y>
                    o.x = x(:)';
                    o.y = y(:)';
                case 3 %<nr x y>
                    o.nr = nr;
                    o.x = x(:)';
                    o.y = y(:)';
                case 4,
                    o.nr = nr; o.id=o.nr;
                    o.x  = x(:)';
                    o.y  = y(:)';
                    o.z  = z(:)';
                otherwise
                    error('%s constructor, wrong number of arguments(%d)',nargin);
            end        
        end
        function o = shift(o,dx,dy)
            o.x=o.x+dx;
            o.y=o.y+dy;
        end
        function o = shiftTo(o,xn,yn)
            o.x = xn;
            o.y = yn;
            o = hydObj(o.nr,o.x,o.y,o.z);
        end
        function xc  = get.xc(o), xc=mean(o.x); end
        function yc  = get.yc(o), yc=mean(o.y); end
        function box = get.box(o); box=[min(o.x),min(o.y),max(o.x),max(o.y)]; end
        function edge= get.edge(o)
            N = length(o.x);
            edge(N).x = [o.x(end) o.x(1)];
            edge(N).y = [o.y(end) o.y(1)];
            for i=(N-1):-1:1
                edge(i).x = [o.x(i) o.x(i+1)];
                edge(i).y = [o.y(i) o.y(i+1)];
                edge(i).L = sqrt(diff(edge(i).x).^2 + diff(edge(i).y).^2);
            end
        end
        function perimeter = get.perimeter(o)
            perimeter= sum([o.edge.L]);
        end

        function area=get.area(o)
            dx=o.x(:)-o.xc; dy=o.y(:)-o.yc;
            area=sum(cross([dx             dy           zeros(size(dx))],...
                           [dx([2:end 1]) dy([2:end 1]) zeros(size(dx))]))/2;
            area=abs(area(end));
        end
        function o=rotate(o,xc,yc,alfa)
            % rotate object around point xc,yc over alfa counter clock wise
            [o.x,o.y]=mf_rotate(o.x-xc,o.y-yc,0,0,alfa);
            o = hydObj(o.nr,o.x,o.y,o.z);
        end
        function plot(o,varargin), plot(o.x,o.y,cl,varargin{:}); end
        function whoami(o); fprintf('I am a <<%s>>\n',class(o)); end
    end
end

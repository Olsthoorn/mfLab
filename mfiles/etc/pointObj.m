classdef pointObj
    % svd tests according to www.uwlax.edu/faculty/will/svd/
    properties
    x  = []
    y  = []
    z  = []
    c  = '';
    clrs = 'brgkmcy';
    end
    properties (Dependent=true)
        P
    end
    methods
        function o = pointObj(x,y,z)
            if nargin==0; return; end
            
            o.x = x;
            o.y = y;
            o.z = z;
            
            o.y(end+1:numel(o.x))=o.y(end); o.y(length(o.x)+1:end)=[];
            o.z(end+1:numel(o.x))=o.z(end); o.z(length(o.x)+1:end)=[];
            
            for it = numel(o.x):-1:1
                o(it).c = mf_color(it,o.clrs);
            end
        end
        function o = circle(o,r,n)
            theta = 2*pi*(0:n-1)/n;
            o.x = r*cos(theta);
            o.y = r*sin(theta);
            o.z = zeros(size(o.x));
            for i=n:-1:1
                o.c(i) = mf_color(i,o.clrs);
            end
        end
        function o = hit(o,A)
            A(:,end+1:min(3,size(A,2))) = [];
            A(end+1:min(3,size(A,1)),:) = [];
            A(:,end+1:3)=0;
            A(end+1:3,:)=0;
            
            for ip = 1:numel(o)
                for i=1:numel(o.x)
                    v = A * [o(ip).x(i); o(ip).y(i); o(ip).z(i)];
                    o(ip).x(i) = v(1);
                    o(ip).y(i) = v(2);
                    o(ip).z(i) = v(3);
                end
            end
        end
        function o = plot(o)
            hold on;
            for ip=1:numel(o)
                plot3([o(ip).x o(ip).x(1)],[o(ip).y o(ip).y(1)],[o(ip).z o(ip).z(1)]);
                for i=1:numel(o.x)
                    plot3(o(ip).x(i),o(ip).y(i),o(ip).z(i),['o' o(ip).c(i)]);
                end
            end
        end
        function o = connect(o,p)
            if numel(o.x) ~= numel(p.x)
                return
            end
            for i=1:numel(o.x)
                plot3([p.x(i),o.x(i)],...
                      [p.y(i),o.y(i)],...
                      [p.z(i),o.z(i)],...
                       o.c(i));
            end
        end
    end
end
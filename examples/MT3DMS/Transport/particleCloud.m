classdef particleCloud
%Folowing an REV (CT4420tranport)
    properties
        NP = NaN;
        dt = 1;
        x  = NaN;
        y  = NaN;
        dx = NaN;
        dy = NaN;
        ux = NaN; % uniform flow 
        uy = NaN; % uniform flow
        a  = NaN; % flow accelration in space
        d  = NaN; %brownian movement step
        hdl= NaN;
        clr= 'b';
    end
    methods
        function PC=particleCloud(x,y,clr,field)
            PC.x=x;
            PC.y=y;
            PC.NP=length(y);
            PC.clr=clr;
            PC.d = field.scale;
            PC.a     = field.a;
            PC.ux    = field.ux;
            PC.uy    = field.uy;
        end
        function PC=plot(PC)
            PC.hdl=plot(PC.x,PC.y,[PC.clr '.'],'markersize',6);
        end
        function PC=show(PC)
            set(PC.hdl,'xdata',PC.x,'ydata',PC.y);
            drawnow;
        end
        function PC=draw(PC)
            PC.hdl=plot(PC.x,PC.y,PC.clr,'linewidth',2);
        end
        function PC=brown(PC) %random step for all object points
            PC.x=PC.x+PC.d*(rand(size(PC.x))-0.5);
            PC.y=PC.y+PC.d*(rand(size(PC.x))-0.5);
            PC=PC.show;
        end
        function PC=flow(PC) % velocity field
            PC.x=PC.x + PC.dt*( + PC.a*PC.x+PC.ux);
            PC.y=PC.y + PC.dt* (- PC.a*PC.y+PC.uy);
            PC=PC.show;
        end
    end
end
    
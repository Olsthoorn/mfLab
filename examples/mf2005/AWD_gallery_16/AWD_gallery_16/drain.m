classdef drain < hydobj
    %DRAIN < hydobj
    %   Detailed explanation goes here
    
    properties
        z
        xshaft
        yshaft
        zshaft
        z_bot
        rw
        Rw
        c0
        cL
        n
        dphi
        PDr
        L
        zone
        I
        ix
        iy
        iz
        LRC
    end
    
    methods
        function obj = drain(name,x,y,z_top,z_bot)
            obj=obj@hydobj(x,y);
            obj.name=name;
            obj.z_top=z_top;
            obj.z_bot=z_bot;
        end
        function obj = rotate(obj,xc,yc,alfa)
            obj=rotate@hydobj(obj,xc,yc,alfa);
            [obj.xshaft,obj.yshaft]=rotate(obj.xshaft-xc,obj.yshaft-yc,alfa);
        end
    end
    
end


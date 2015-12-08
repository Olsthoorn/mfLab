classdef piezometer < hydobj
    %PIEZOMETER < hydobj
    %   Detailed explanation goes here
    
    properties
        z_top
        z_bot
        t
        head
    end
    
    methods
        function obj = piezometer(name,x,y,z_top,z_bot)
            if nargin==0
                name='piezometer';
                x=0;
                y=0;
                z_top=NaN;
                z_bot=NaN;
            end
            obj=obj@hydobj(x,y);
            obj.name=name;
            obj.z_top=z_top;
            obj.z_bot=z_bot;
        end
    end
    
end


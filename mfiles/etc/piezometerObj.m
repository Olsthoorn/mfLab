classdef piezometerObj <hydObj
%PIEZOMETEROBJ class definition of piezometer Objects
%
% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later
    %PIEZOMETER < hydobj
    %   Detailed explanation goes here
    properties
        z_top
        z_bot
        t
        h
    end
    methods
        function obj = piezometerObj(name,x,y,z_top,z_bot)
            %PIEZOMETEROBJ constructor for a piezometer Object.
            
            obj.type='piezometer';
            switch nargin
                case 0, return;
                case 2, % interprete piezometer(pieznames(m,1),list(m,[x,y,z_top,z_bot])
                    data=x; % interprete second argument as list [x,y,z_top,z_bot]
                    m=size(data,1);
                    if m~=numel(name),
                        error('%s constructor: Nr of names(%d) ~= length of of props(%d)',class(obj),numel(name),m);
                    end
                    for i=m:-1:1
                        obj(i)=piezometerObj(name{i},...
                            data(i,1),data(i,2),data(i,3),data(i,4));
                    end
                case 5                   
                    obj.name=name;
                    obj.x=x;
                    obj.y=y;
                    obj.z_top=z_top;
                    obj.z_bot=z_bot;
                otherwise
                    error('%s constructor: Wrong number of arguments.',class(obj));
            end
        end
        function obj=seth(obj,t,h)
            if all(t<datenum(1899,12,30)), t=t+datenum(1899,12,30); end % convert from Excel
            obj.t=t; 
            obj.h=h;
        end
        function obj=chart(obj,clrlt,t1,t2)
            hold on
            if nargin==1
                clrlt='ob-';
            end
            if nargin==4
                plot(obj.t(obj.t>=t1 & obj.t<=t2),obj.h(obj.t>=t1 & obj.t<=t2),clrlt);
            else
                plot(obj.t,obj.h,clrlt);
            end
        end
        function obj=plot(obj,cl)
                plot(obj.x,obj.y,cl);
        end
    end
end


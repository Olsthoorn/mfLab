classdef ADVmethod
%ADVMETHOD class definition for advection methods
%
% TO 123112
%
% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

    properties
        value
        name
    end
    methods
        function o=ADVmethod(m)
           %ADVMETHOD --- constructor for ADVmethod objects
           % TO 121231
           
            switch(m)
                case {-1,'TVD'}, o.value=-1; o.name='Ultimate, TVD';
                case { 0,'FDM'}, o.value= 0; o.name='FDM, standard';
                case { 1,'MOC'}, o.value= 1; o.name='MOC';
                case { 2,'MMOC'},o.value= 2; o.name='MMOC';
                case { 3,'HMOC'},o.value= 3; o.name='HMOC';
                otherwise
                    o.value= -1;
                    o.name = 'FDM, standard';
            end
        end
    end
end

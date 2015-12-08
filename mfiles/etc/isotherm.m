classdef isotherm
%ISOTHERM determines sorption isotherm (class def)
%
% USAGE:
%    none

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later
    properties
        value
        name
    end
    methods
        function o=isotherm(m)
            %ISOTHERM isotherm object constructor
            %
            % o=isotherm(m) --- m if sorption type
            %
            % TO 130428
            
            switch(m)
                case {0,'none'      }, o.value=-1; o.name='no sorption';
                case {1,'linear'    }, o.value= 0; o.name='Linear';
                case {2,'freundlich'}, o.value= 1; o.name='Freundlich';
                case {3,'langmuir'  }, o.value= 2; o.name='Landmuir';
                case {4,'kinetic'   }, o.value= 3; o.name='First-order kinetic';
                case {5,'ddmt'      }, o.value =4; o.name='Dual-domain mass transfer';
                case {6,'ddmts'     }, o.value =5; o.name='Dual-domain mass transfer + sorption';
                otherwise
                    o.value= 0;
                    o.name = 'none';
            end
        end
    end
end

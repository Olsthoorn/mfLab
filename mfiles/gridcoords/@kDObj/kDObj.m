classdef kDObj
    % kDObj general hydrological object
    % object encountered in the WEL file of the Dutch Natioal Hydrologic
    % Instrument (www.NHI.nu), specifying transmissivity kD for a x,,y
    % location and model-layer depth.
    % It is loaded and can then be questioned like with any struct
    % TO 120425 120815
    properties
        % Properties made as similar to wellObj as possible, many may never be
        % necessary but could eventually be usefull.
        nr,        % nr may be same as id (mostly not necessary)
        id,        % may be same as id (mostly not necessary)
        name,      % any practical name of this kD
        longName,  % any long name, may be the pumping test location
        x, y,      % coordianates
        z,         % z is top and bottom of aquifer, zTop gr
        ztop,      % zTop surface elevation
        DZ,        % layer thickness (for transmissivity)
        ix, iy, iLay, idx, LRC,
        parent,     % perhaps and aquifer name or so
        remark =[],
        created,
        code,      % for whatever purpose
        value,     % transmissivity value
        type = 'kDObj',
        UserData   % for the user to add any data
    end
    methods
        function obj=kDObj()
            switch nargin
                case 0,
                    o.created = now;
                otherwise
                    error('%s constructor, wrong number of arguments(%d)',nargin);
            end        
        end
                function X   = X(o)
            % get first x coordinate of wells
            for io=numel(o):-1:1
                X(io) = o(io).x(1);
            end
        end
        function Y = Y(o)
            % get first y coordinate of wells
            for io=numel(o):-1:1
                Y(io) = o(io).y(1);
            end
        end
        function Z_top = Z_top(o)
            % get highest z-coordinate of well
            for io=numel(o):-1:1
                Z_top(io)= o(io).z(1);
            end
        end
        function Z_bot = Z_bot(o)
            % get lowest z-coordinate of well
            for io=numel(o):-1:1
                Z_bot(io)=o(io).z(end);
            end
        end
        function Idx = Idx(o)
            % Get first index of every wellObj in a vector
            % USAGE: Idx= wellObj.Idx;
            for io=numel(o):-1:1
                Idx(io)=o(io).idx(1);
            end
        end
        
        function Ix= Ix(o)
            for io=numel(o):-1:1
                Ix(io) = o(io).ix(1);
            end
        end
        
        function Iy= Iy(o)
            for io=numel(o):-1:1
                Iy(io) = o(io).iy(1);
            end
        end

    end
end

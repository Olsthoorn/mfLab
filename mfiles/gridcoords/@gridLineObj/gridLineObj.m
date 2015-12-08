classdef gridLineObj
    %GRIDLINEOBJ class def to represent a line that intersets a finite difference grid.
    % It takes the following forms, where the number of points is undefined. That is
    % the line may consist of any number of points. It may have points that  are 
    % somewhere in a cell, like the begin and end points or intermediate
    % points in case the original line  was polyline, and it has all points
    % where the line intersets a cell face. Any line piece therefore, has
    % at least two coordinates and is completely inside a cell. If urther has
    % the cell indices in the grid, including the global index in the grid.
    % It also has its length
    % The objec't user data can be used to store additional information.
     properties
         x; y; z; zR % coordinates of line piece, all two values
                     % zR is realative z-coordinates
         ix; iy; iz; % cell indices
         L           % total length of line within cell
         idx         % cell global index
         iLay        % layer number if model cell layer else 0
         iCbd        % aquitard number if aquitard else 0
         isLay       % true if isLayer
         isCbd       % true if isCBD
     end
     properties (Dependent=true)
         xm; ym; zm; zRm % center of line piece
     end
     methods
        function o = gridLineObj(gr,pline)
            %GRIDLINEOBJ constructor for line through grid
            %
            % Usage:
            %   gridLineObj= gridLineObj(gr,pline) % generate a gridLineObject
            
            if nargin<2; return; end
            
            % Make sure pline has at least 2 points
            if size(pline,1)<2, error('%s: %s must have at least to points',mfilename,class(o)); end

            % 
            % Make sure pline is a 3D line
            if size(pline,2)~=3, error('%s: %s must be 3D [x y z]',mfilename,class(o)); end
            
            % Add column ot hold relative z coordinates
            pline = [pline pline(:,end)];
            
            % Get array of line piece objects from cutting pline through
            % the grid. Each line piece is inside a model cell (no
            % aquitard). It has all properties.
            o=gr.lineObjects(pline);
        end
        function  o = mergeDoubles(o)
            %MERGEDOULBLE -- mergens the gridLineObj in the same cell.
            % USAGE: P = P.mergeDoulbe();
            %   P is an array of gridLineObj. gridLineObj targeting the
            %   same cell will be merged. That is, their x, y and z
            %   coordinates are concatenated in the order of the P-elements
            %   in the P array.
            for io=numel(o)-1:-1:1
                if o(io).idx == o(io+1).idx
                    o(io).x = [o(io).x; o(io+1).x];
                    o(io).y = [o(io).y; o(io+1).y];
                    o(io).z = [o(io).z; o(io+1).z];
                    o(io+1) = [];
                end
            end
            
        end
        function xm = get.xm(o), xm = mean(o.x(:)); end
        function ym = get.ym(o), ym = mean(o.y(:)); end
        function zm = get.zm(o), zm = mean(o.z(:)); end
    end
end
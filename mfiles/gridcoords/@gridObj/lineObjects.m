function [grLines,L]=lineObjects(o,pline,zRel) % a method of grid
%GRIDOBJ/LINEOBJECTS --- generates a vector of gridLineObj that represent the
% intersectinon of the line in pline with the grid.
% Each gridLineObj contains the info of the piece of the total line
% that intersect a cell. It has the cell number, the entry and exist points
% and the length of the line inside the intersected cell.
% This method of the gridObj is generally called by te constructor of
% gridLineObj.
%
% USAGE:
%     [grLineBuf,L] = grid.lineObjects(pline3D);
%     [grLineBuf,L] = grid.lineObjects(pline2D,zRel);
%
%  lineObjects returns grLineBuf, which is an array with lineObjects that
%  contain the information of each model cell intersected by the pline;
%
% L is cell array size pline holding the idx of the cells cut by each line
% section (not vertex).
%
%  pline3D is a vector with 3 columns, corresponding to x, y and z of the line
%  pline2D is a vector with 2 columns, corresponding to x and y of the line
%  zRel is only used in conjuncties with pline2D. It is the relative
%  z-coordinate of the pline in the grid. Relative coordinate a.b means
%  layer a, fraction b. I.e. 3.5 is the center of layer 4, counted from the
%  top of the model.
%
% TO 130826 complete overhaul now using abolute z to allow straight
% khettaras to be defined or straight drain betwween the vertices
% irrespective of the intemediate grid elevation.
% You can use relative coordinates for all vertices though. But between
% them absolute z-coordinates are used as derived for the cells in which the
% vertices happen to be (old version saved as lineObjects2 (Not in SVN).
%
% TO 121119 130626 130826

% vertices 

    % preallocate NBUF grLineBuf objects to have a buffer
    NBUF = 100; ibuf = 0; grLineBuf(NBUF,1) = gridLineObj();

    % remove parts of line totally outside the grid, this concerns the
    % points who and their both neighbours are all outside the grid

    % to know which points are outside the model
    IDX = xyzindex(pline(:,[1 2]),o);
    
%     if any(isnan(idx))
%         J = find(isnan(idx));
%         nanChk = [NaN; idx; NaN];
%         for j=numel(J):-1:1
%             i=J(j);
%             if all(nanChk(i+(0:2))), pline(i,:)=[]; idx(i)=[]; end
%         end
% 
%         if size(pline,1)<2
%             grLines = grLineBuf(end); % this is an empty gridLineObj
%             return;
%         end
%     end
    
    isZrelative = nargin>2;
    
    if isZrelative
        % relative z-coordinate(s) provided, then change the grid o to
        % a grid using relative z relative coordinates
        
        gr = o;         % remember the current absolute grid
        o  = gridObj(gr.xGr,gr.yGr,-(0:gr.Nz),gr.LAYCBD,gr.MINDZ,gr.AXIAL);
        
        zRel(length(zRel):size(pline,1))= zRel(end); zRel=zRel(1:size(pline,1));
        zRel  = min(o.Nlay,zRel);
        pline = [pline(:,1:2) -zRel(:)];
        if any(zRel<0 | zRel>o.Nlay)
            error(['Relative zCoordiate(s), zRel must be between\n',...
                   'base and top of system).\n',...
                   'zRel = (iLay-1) + fraction of aquiferThickness']);
        end
    end

    % continue as if all z coordinates are absolute
    
    if size(pline,1)<2
        error('%s: pline needs at least two vertices',mfilename);
    end
    
    % length of line pieces
    dxyz = diff(pline(:,[1 2 3]),1,1);

    nGrLinesPrev = 1; % no gridLines yet. To remember cells pertaining to each line section
    
    % Don't use relative coordinates, but follow the line straight through the grid
    for iL = 1:size(pline,1)-1 % includes first and last point
        dx=dxyz(iL,1);
        dy=dxyz(iL,2);
        dz=dxyz(iL,3);

        % first find the column hit by the line in horizontal projection
        % i.e. don't assume the vertex is in the grid at all, yet the line
        % may intersect the grid further down.
        lambdaX= (o.xGr(:)' - pline(iL,1))/dx;
        lambdaY= (o.yGr(:)' - pline(iL,2))/dy;
                
        % combine all lambdas,
        %lambdaXY = roundn([lambdaX lambdaY],6);
        lambdaXY = [lambdaX lambdaY];
        
        % make sure to inlcude the first and the last point by using lambdaXY 0 and 1
        % This is always necessary as most lines start and end somewhere
        % inside a cell and almost never at s cell boundary.
        % if IDX(iL  ) ~= NaN, the point is inside the model,  so add 0 to lambda, likewise
        % if IDX(iL+1) ~= NaN, the end point is inside the model, add 1 tolambda
        if ~isnan(IDX(iL  )), lambdaXY = [0 lambdaXY]; end %#ok start point inside model
        if ~isnan(IDX(iL+1)), lambdaXY = [lambdaXY 1]; end %#ok end   point insize model
        
        % round to prevent missing hits in numerical comparisons
        lambdaXY = unique(lambdaXY(lambdaXY>=0 & lambdaXY<=1));
        
        if isempty(lambdaXY), continue; end
        
        % these lambdas refer the hits of this line piece with xGr and yGr
        % this leaves at least [0 1], which refer to points inside the model
        % but not on on of the cell faces, i.e. one of the vertical cell faces.
        % Juse [0 1] is a line totally within a single cell (projected to
        % the xy plane, not considering the z.
        % Hence, lambdaXY can never be empty and it should even have at least
        % two points
        if length(lambdaXY)<2
            continue;
        end
        
        % Mids of the line pieces between intersections with xGr and yGr
        % guarantee points to be within a cell (projected on the xy plane).
        lambdaXYc=0.5*(lambdaXY(1:end-1)+lambdaXY(2:end));

        % the mids of the cell intersections in the XY plane.
        % xc and yc should be at least of length 1.
        xc=pline(iL,1)+lambdaXYc*dx;
        yc=pline(iL,2)+lambdaXYc*dy;        

        % we continue with cleand-up lambda and 
        % x,y and z where line intersects a cell face
        xFace    = pline(iL,1)+lambdaXY*dx; 
        yFace    = pline(iL,2)+lambdaXY*dy;
        z_xyFace = pline(iL,3)+lambdaXY*dz; % z at xFace,yFace not at zFace
              
        % z_xyFace is the z coordinate at the vertical faces intersected by
        % the original line section. To find which cells it intersects we
        % have to look within every project xy cell which stories in the
        % vertical it intersects, because these stories (z-elevations)
        % differ between cells unless all layers are uniform.
        dxFace= diff(xFace);
        dyFace= diff(yFace);
        dzFace= diff(z_xyFace);

        % indices of each intersected horizontally projected cels:
        [ixc,iyc] = xyzindex(xc,yc,o);

        % the task now is to generate grLineObjects for each intersected
        % cell in the xy plane given by ixc and iyc. That is, within the
        % column that belongs to this horizontally projected cell.
        % This is doen by intersection of the line with the horizontal planes
        % at the cell top and bottom elevations within this column.
        % ixc,iyc should not be NaN as points outside the horizotnal projection
        % of the model have already be eliminated.
        
        % For each horizontally projected cell intercectec by our line section:
        for ic = 1:numel(ixc)
            
            if isnan(ixc(ic)) || isnan(iyc(ic))
                continue;
            end
            
            % cutting line from first hit vertically
            zGr = o.Z(iyc(ic),ixc(ic),:);
            zGr = zGr(:)';
            
            % if both points are below or above the cells in the column,
            % continume, ignore this portion of the line
            if all(z_xyFace(ic+[0 1])<zGr(end)) || all(z_xyFace(ic+[0 1])>zGr(1))
                continue;
            end
            
            % IDZ is NaN for points outside the grid (within this column)
            IDZ = xyzindex(z_xyFace(ic+[0 1]),zGr);
            
            lambdaZ = (zGr-z_xyFace(ic))/dzFace(ic);
            
            % We add zero and 1 because the start and end z could fall
            % totally with a single cell if IDZ~=NaN
            
            if ~isnan(IDZ(  1)), lambdaZ = [0 lambdaZ]; end %#ok, point inside the grid
            if ~isnan(IDZ(end)), lambdaZ = [lambdaZ 1]; end %#ok, point inside the grid
            lambdaZ = unique(lambdaZ);
            lambdaZ = lambdaZ(lambdaZ>=0 & lambdaZ<=1);
            if numel(lambdaZ)<2, continue; end
            lambdaZc= 0.5*(lambdaZ(1:end-1)+lambdaZ(2:end));
            
            % lambdaZ cannot be empty as it as at least [0 1] and so
            % lambdaZc has at least one point referring to a cell.
            
            zc = z_xyFace(ic) + lambdaZc * dzFace(ic); % at least 1 point [0.5]            
            iz = xyzindex(zc(:),zGr(:));            
            % however, the first and or the last zc may be outside the model
            
            if isnan(iz(  1))
                lambdaZ(  1) = []; iz(  1)=[];
                if isempty(iz), continue; end
            end
            if isnan(iz(end))
                lambdaZ(end) = []; iz(end)=[];
                if isempty(iz), continue; end
            end
            
            if any(isnan(iz))
                error('iz may not contain NaN''s at this stage, check logic');
            end
            
            % iz is not empty at this point, so determine the grLines
            for i=1:numel(iz)
                if ibuf==NBUF
                   ibuf = 1;
                   if ~exist('grLines','var') || isempty(grLines)
                       grLines = grLineBuf;
                   else
                       grLines = [grLines; grLineBuf]; %#ok
                   end
                    clear('grLineBuf');
                    grLineBuf(NBUF,1) = gridLineObj();
                else
                    ibuf=ibuf+1;
                end                

                grLineBuf(ibuf).x     = xFace(ic)+dxFace(ic)*lambdaZ(i+[0 1]);
                grLineBuf(ibuf).ix    = ixc(ic);

                grLineBuf(ibuf).y     = yFace(ic)+dyFace(ic)*lambdaZ(i+[0 1]);
                grLineBuf(ibuf).iy    = iyc(ic);

                z                     = z_xyFace(ic)+dzFace(ic)*lambdaZ(i+[0 1]);     % absolute coordinates
                grLineBuf(ibuf).z     = z;
                grLineBuf(ibuf).iz    = iz(i);

                Idx                   = cellIndex(ixc(ic),iyc(ic),iz(i),o);
                grLineBuf(ibuf).idx   = Idx;

                grLineBuf(ibuf).zR    = roundn((o.Z(Idx)-z)/o.DZ(Idx)+iz(i)-1,4);
                
                % There may be some issues here with LAYCBD (if confining
                % beds are present in the grid).
                % We have tried to stay away from LAYCBD by using gr.Z, so
                % that iz is referring to both LAY and CBD layers.
                % We add these lines to indicate whether iz is a layer or a
                % confining bed.
                grLineBuf(ibuf).isLay =  o.isLay(iz(i));
                grLineBuf(ibuf).isCbd = ~o.isLay(iz(i));
                grLineBuf(ibuf).iLay  = sum( o.isLay(1:iz(i))) *  o.isLay(iz(1));
                grLineBuf(ibuf).iCbd  = sum(~o.isLay(1:iz(i))) * ~o.isLay(iz(i));
                grLineBuf(ibuf).L     = sum(sqrt(diff(grLineBuf(ibuf).x,1,2).^2+...
                                                 diff(grLineBuf(ibuf).y,1,2).^2+...
                                                 diff(grLineBuf(ibuf).z,1,2).^2));
                % merge grLines if idx is the same in consecutive ones
                if ibuf>1
                    if Idx==grLineBuf(ibuf-1).idx
                        grLineBuf(ibuf-1).x(end) = grLineBuf(ibuf).x(end); 
                        grLineBuf(ibuf-1).y(end) = grLineBuf(ibuf).y(end);
                        grLineBuf(ibuf-1).z(end) = grLineBuf(ibuf).z(end);
                        grLineBuf(ibuf-1).zR(end)= grLineBuf(ibuf).zR(end);
                        grLineBuf(ibuf-1).L      = grLineBuf(ibuf-1).L + grLineBuf(ibuf).L;
                        ibuf = ibuf-1;
                    end
                elseif exist('grLines','var') && numel(grLines)>1
                    if Idx==grLines(end).idx
                        grLines(end).x(end) = grLineBuf(ibuf).x(end); 
                        grLines(end).y(end) = grLineBuf(ibuf).y(end);
                        grLines(end).z(end) = grLineBuf(ibuf).z(end);
                        grLines(end).zR(end)= grLineBuf(ibuf).zR(end);
                        grLines(end).L      = grLines(end).L + grLineBuf(ibuf).L;
                        ibuf = ibuf-1;
                    end
                else
                    continue;
                end
            end
        end
        % how many gridLines?
        if ~exist('grLines','var')
            nGrLines=0;
        else
            nGrLines = numel(grLines);
        end
        if nGrLinesPrev<nGrLines % then some in buffer
            Ig = nGrLinesPrev:nGrLines; % so many direct
            Ib = 1:ibuf;                       % these from buffer
        else
            Ig = [];                           % all in buffer
            Ib = max(nGrLinesPrev-nGrLines,1):ibuf;
        end
        if isempty(Ig)
            L{iL} = [grLineBuf(Ib).idx]; %#ok
        else
            L{iL} = [[grLines(Ig).idx] [grLineBuf(Ib).idx]]; %#ok
        end
        nGrLinesPrev = nGrLines+ibuf;
    end
    
    if ~exist('grLines','var')
        grLines = grLineBuf(1:ibuf);
    else
        grLines = [grLines; grLineBuf(1:ibuf)];
    end
    
    if ~isempty(grLines)
        grLines([grLines.isCbd]) = [];
    end
    
    % if z is relative compute the absolute coordinates
    % else compute the relative coordinates
    if isZrelative
        for i=numel(grLines):-1:1
          zR = -grLines(i).z;    frac_zR = zR -grLines(i).iz+1;
          idx=  grLines(i).idx;
          
          grLines(i).z  =  gr.Z(idx)-gr.DZ(idx)*frac_zR;
        end
    else % compute the relative coordinates
        for i=numel(grLines):-1:1
            z  = grLines(i).z;
            idx= grLines(i).idx;
            
            grLines(i).zR = grLines(i).iz-1+(o.Z(idx)-z)/o.DZ(idx);
        end
    end
end

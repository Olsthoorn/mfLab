function varargout=xyzindex(varargin)
%XYZINDEX computes the cell indices of points in 1D, 2D or 3D grid
%    and optionally the relative coordinates as well
%    ix,iy,iz are cell indices, xR,yR and zR are relative coordinates inside
%    the cells.
%    The number of input requesting coordinates, xW,yW and zW should be the
%    same size. If not they are made so by extendign the shortest vectors.
%    You may also combine these inputs to array of 2 or 3 columns in 2 and
%    3D. Notice that the length of xGr, yGr and zGr are generally
%    different. Therefore, they cannot be combined into arrays with 2 or 3
%    columns and must be given as vectors, or alternatively as a gridObj.
%
% EXAMPLES: (all possible input combinations)
%  without gridObj:
%   [ix,xR,xW]                   = xyzindex( xW       ,xGr);
%   [ix,iy,xR,yR,xW,yW]          = xyzindex( xW,yW    ,xGr,yGr);
%   [ix,iy,xR,yR,xW,yW]          = xyzindex([xW yW]   ,xGr,yGr);
%   [ix,iy,iz,xR,yR,zR,xW,yW,zW] = xyzindex( xW,yW,zW ,xGr,yGr,zGr);
%   [ix,iy,iz,xR,yR,zR,xW,yW,zW] = xyzindex([xW,yW,zW],xGr,yGr,zGr);
%  with one output argument:
%    ix                          = xyzindex( xW       ,xGr);
%    idx                         = xyzindex( xW,yW    ,xGr,yGr);
%    idx                         = xyzindex([xW yW]   ,xGr,yGr);
%    Idx                         = xyzindex( xW,yW,zW ,xGr,yGr,zGr);
%    Idx                          = xyzindex([xW,yW,zW],xGr,yGr,zGr);
%  with gridObj:
%   [ix,xR,xW]                   = xyzindex( xW       ,gr);
%   [ix,iy,xR,yR,xW,yW]          = xyzindex( xW,yW    ,gr);
%   [ix,iy,xR,yR,xW,yW]          = xyzindex([xW yW]   ,gr);
%   [ix,iy,iz,xR,yR,zR,xW,yW,zW] = xyzindex( xW,yW,zW ,gr);
%   [ix,iy,iz,xR,yR,zR,xW,yW,zW] = xyzindex([xW,yW,zW],gr);
%  with one output argument
%    ix                          = xyzindex( xW       ,gr);
%    idx                         = xyzindex( xW,yW    ,gr);
%    idx                         = xyzindex([xW yW]   ,gr);
%    Idx                         = xyzindex( xW,yW,zW ,gr);
%    Idx                         = xyzindex([xW,yW,zW],gr);
%
% ix,iy,iz are cell indices
% xR,yR,zR are relative coordinates of xW, where xR=5.3 is xGr(6)+0.3dx(5),
%          i.e. at 0.3 of cell 5 relative coordinates start at xR=0 (xW==xGr(1)
%          to xR==Nx, i.e. xW==xGr(end)
% xGr,yGr,zGr are grid coordinates, Z may be a complete Z grid of size [Nx,Ny,Nz+1]
% xW,yW,zW are world coordinates
%
% With one output argument the global indices Idx are produced in the 3D case.
% In the 2D case, idx is the global index in layer one. That is, within the
% layer. To get it for layer iLay use Idx = idx + Nx*Ny*(iLay-1).
%
% with argument gr as gridObj, LAYCBD (confining beds) is taken into
% account. LAYCBD is a property of the gridObj (see gr.LAYCBD)
% in CBD layers, the relative coordinate runs from -1 to 1 while the
% model layer number is the same as that of the overlaying layer.
%
% see also testxyzindex cellIndex cellIndices linegrid inpolygon inpoly IdxMatlab2Modflow
%
% TO 090317 090509 100510 110818 120513 121128 130127
% TO 130826 update, new help [xGr,yGr] and [xGr,yGr,zGr] now illegal, must
%           now be vectors.

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

errMsg =['without a grid, the last input argument must be a vector.\n',...
         'If there is no grid, the last argument must be a vector because\n',...
         'xGr,yGr and zGr are generally of different size and cannot be\n',...
         'specified as an array with more than one column.'];

[gr,varargin] = getType(varargin,'gridObj',[]);

if isempty(varargin)
    error('No inputs.');
end

if isempty(gr)
    firstIsMultiCol = size(varargin{1},1)==1 && ( ...
            size(varargin{1},2)==2 && numel(varargin)==3 || ...
            size(varargin{1},2)==3 && numel(varargin)==4 ...
        );
else
    firstIsMultiCol = size(varargin{1},1)==1 && size(varargin{1},2)<=3;
end

for i=1:numel(varargin)
    varargin{i}=varargin{i}(:); % assert vectors are all column vectors
end

if firstIsMultiCol
    varargin{1} = varargin{1}';
end

if isempty(gr)
    if numel(varargin)<2, error(errMsg); end
else
    if numel(varargin)<1, error(errMsg); end
end

if ~isempty(gr)  % varargin contains a grid
    Nx = gr.Nx;
    Ny = gr.Ny;    
    switch numel(varargin)  % number of non-grid inputs
        case 1
            % argument is an array with 1, 2 or 3 vertical vectors
            switch size(varargin{1},2)
                case 1  % [ix,xR,xW] = xyzIndex(xW,gr)
                    xW = varargin{1}(:,1); [varargout{1},varargout{2}] = getindices(xW,gr.xGr);
                    if nargout>2
                        varargout(3) = {xW};
                    end
                case 2  % [ix,iy,xR,yR,xW,yW] = xyzIndex([xW yW],gr)
                    xW = varargin{1}(:,1); [varargout{1},varargout{3}] = getindices(xW,gr.xGr);
                    yW = varargin{1}(:,2); [varargout{2},varargout{4}] = getindices(yW,gr.yGr);
                    if nargout>4
                        varargout(5:6) = {xW, yW};
                    end
                case 3  % [ix,iy,iz,xR,yR,zR,xW,yW,zW] = xyzIndex([xW yW zW],gr)                    
                    xW = varargin{1}(:,1); [varargout{1},varargout{4}] = getindices(xW,gr.xGr);
                    yW = varargin{1}(:,2); [varargout{2},varargout{5}] = getindices(yW,gr.yGr);
                    zW = varargin{1}(:,3); [varargout{3},varargout{6}] = zIndices(zW,varargout{1},varargout{2},gr);
                    if any(gr.LAYCBD>0)
                        [varargout{3},varargout{6}] = zlay2model(varargout{3},varargout{6},gr);
                    end
                    if nargout>6
                        varargout(7:9) = {xW, yW, zW};
                    end
                otherwise
                    error('%s: don''t understand combination of inputs',mfilename);
            end
        case 2 % [ix,iy,xR,yR,xW,yW] = xyzIndex(xW,yW,gr)
            xW = varargin{1}; [varargout{1},varargout{3}] = getindices(xW,gr.xGr);
            yW = varargin{2}; [varargout{2},varargout{4}] = getindices(yW,gr.yGr);
            if nargout>4
                varargout(5:6) = {xW, yW};
            end
        case 3  % [ix,iy,iz,xR,yR,zR,xW,yW,zW] = xyzIndex(xW yW zW,gr)                    
            xW = varargin{1}; [varargout{1},varargout{4}] = getindices(xW,gr.xGr);
            yW = varargin{2}; [varargout{2},varargout{5}] = getindices(yW,gr.yGr);
            zW = varargin{3}; [varargout{3},varargout{6}] = zIndices(zW,varargout{1},varargout{2},gr);
            if any(gr.LAYCBD>0)
                [varargout{3},varargout{6}] = zlay2model(varargout{3},varargout{6},gr);
            end
            if nargout>6
                varargout(7:9) = {xW, yW, zW};
            end
        otherwise
            error('%s: don''t understand combination of inputs',mfilename);
    end
    
else % arguements do not contain a grid object
    switch numel(varargin) % number of numeric arguments
        case 2, % xyzIndex(xW,xGr)
            try
                xW = varargin{1}; [varargout{1},varargout{2}] = getindices(xW,varargin{2});
                if nargout>2
                    varargout(3) = {xW};
                end
            catch ME
                fprintf('Invalid input combination, two args varargin must be {xY,xGr} or {xW, gr}\n');
                throw(ME);
            end
        case 3 % xyzIndex([xW yW],xGr,yGr)
            try 
               xW = varargin{1}(:,1); [varargout{1},varargout{3}] = getindices(xW,varargin{2});
               yW = varargin{1}(:,2); [varargout{2},varargout{4}] = getindices(yW,varargin{3});
               if nargout>4
                   varargout(5:6) = {xW, yW};
               end
               Nx = length(varargin{2})-1;
               Ny = length(varargin{3})-1;
            catch ME
                fprintf('%s: Don''t understand input combination\n',mfilename);
                throw(ME);
            end
        case 4 % xyzIndex(xW,yW,xGr,yGr) 
            try
                if ~firstIsMultiCol
                   xW = varargin{1}; [varargout{1},varargout{3}] = getindices(xW,varargin{3});
                   yW = varargin{2}; [varargout{2},varargout{4}] = getindices(yW,varargin{4});
                   if nargout>4
                       varargout(5:6) = {xW, yW};
                   end
                   Nx = length(varargin{3})-1;
                   Ny = length(varargin{4})-1;
                else % xyzIndex([xW yW zW],xGr yGr zGr)
                   xW = varargin{1}(:,1); [varargout{1},varargout{4}] = getindices(xW,varargin{2});
                   yW = varargin{1}(:,2); [varargout{2},varargout{5}] = getindices(yW,varargin{3});
                   zW = varargin{1}(:,3); [varargout{3},varargout{6}] = getindices(zW,varargin{4});
                   if nargout>6
                       varargout(7:9) = {xW, yW, zW};
                   end
                   Nx = length(varargin{2}-1);
                   Ny = length(varargin{3}-1);
                end
            catch ME
                fprintf('%s: Don''t understand input combination\n',mfilename);
                throw(ME);
            end
        case 5
            error('%s: Illegal input combination',mfilename);               
        case 6 % xyzIndex(xW,yW,zW,xGr,yGr,zGr)
            try
                   xW = varargin{1}; [varargout{1},varargout{4}] = getindices(xW,varargin{4});
                   yW = varargin{2}; [varargout{2},varargout{5}] = getindices(yW,varargin{5});
                   zW = varargin{3}; [varargout{3},varargout{6}] = getindices(zW,varargin{6});
                   if nargout>6
                       varargout(7:9) = {xW, yW, zW};
                   end
                   Nx = length(varargin{4}-1);
                   Ny = length(varargin{5}-1);
            catch ME
                fprintf('%s: Illegal input combination.\n',mfilename);
                throw(ME);
            end
        otherwise
            error('%s: Illegal input combination',mfilename);
    end
end

% Make sure the length of the outputs are the same
if exist('varargout','var')
    L   = cellfun(@length,varargout);
    Lmax= max(L);

    for i=length(L):-1:1
        if length(varargout{i})<Lmax
            varargout{i}(end+1:Lmax,1)=varargout{i}(end);
        end
    end

    if nargout==1, % then compute Idx (global index)   
        switch length(varargout)
            case 2
                Idx = varargout{1};
            case 4
                Idx = Ny*(varargout{1}-1) + varargout{2};
            case 6
                Idx = Ny*Nx*(varargout{3}-1) + Ny*(varargout{1}-1) + varargout{2};
        end
        varargout{1} = Idx;
    end
else
    varargout{1} = NaN;
end

end

function [Ix,xR] = getindices(xW,xGr)
% GETINDCES get cell indices for all values of xW
% USAGE: Ix = getindices(xW,xGrid)
    xR = interp1(xGr(:),1:numel(xGr),xW(:));
    Ix = min(numel(xGr)-1,floor(xR));
    Ix(isnan(xR)) = NaN;
    xR = xR-Ix;
end

function [Iz,zR] = zIndices(zW,ix,iy,gr)
    % Get zW indices in case we may have a full 3D zGr array
    for j=numel(ix):-1:1
        [Iz(j,1),zR(j,1)] = getindices(zW(j),gr.Z(iy(j),ix(j),:));
    end
end

function [iz,zR] = zlay2model(iz,zR,gr)
    % [iz,zR] = zlay2model(iz,zR,gr) --- take LAYCBD into account if present.
    % Also assign relative coordinates if there are confining beds.
    % Specifically for MODPATH (Version 6):
    % The relative coordinate within a confining beds is -1 at its bottom and 0 at its top.
    % Its layer index equals that of the layer on top of it, in which the
    % relative coordinate varies between 0 and 1.
    
    if all(gr.LAYCBD==0), return; end
    
    % Convert layers to model layers and confining beds
    IZ = [ones(gr.Nlay,1) gr.LAYCBD]';
    IZ = reshape(cumsum(IZ(:)),[2,gr.Nlay])';
    
    IC = [(1:gr.Nlay)' IZ(:,2)]; IC(gr.LAYCBD==0,:)=[];
    IZ = [(1:gr.Nlay)' IZ(:,1)];
    
    IZ = sortrows([IZ;IC],2); % size is Nz,2
    
    % reduce the relative zW coordinate for cbd layers by 1.0
    for i=1:gr.Ncbd
        I = find(iz == IC(i,2));
        if ~isempty(I)
            zR(I) = zR(I)-1.0;
        end
    end 
    iz = IZ(iz,1); % iz is now iLay
end

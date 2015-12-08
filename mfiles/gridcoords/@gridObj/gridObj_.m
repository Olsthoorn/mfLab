classdef gridObj
    %   gr = gridObj(xGr,yGr[,zGr|Z [,'LAYCBD',LAYCBD[,'MINDZ',MINDZ[,'AXIAL',true|false]]])
    %   generates a grid object with properies and methods through which
    %   this object can provide any required information on the grid.
    %   Essential is the difference between layers (= model cell layers,
    %   which are all the defined variables with lay, Lay or LAY in them)
    %   and confining beds (which are all variables with cbd, Cbd or CBD in them).
    %   Confining beds or CBD are not counted as layers.
    %   They are optionally intermediate between each pair of layers
    %   and are connected to the bottom of the layer with which they are associated.
    %   Confining beds have no thickness, so that no elevation is specified
    %   for them in the Z array or zGr input !
    %   The actual configuration of LAY and CBD is specified 
    %   in vector LAYCBD(Nlay,1) in which a value >0 indicates a model layer
    %   with a confining bed below it. Layers without a confining bed
    %   associated below it have LAYCBD=0.
    %   LAYCBD(end)=0 is guaranteed by mfLab.
    %   Z is the elvation of all tops and bottom of layers (excluding
    %   confining beds, which are implied through LAYCBD.
    %        Nlay+1<=size(Z,3)<2NLay
    %   Nz and Nlay are identical. Ncbd is used to address confining bed elevation.
    %   If you require confining beds to have a thickness, then define them
    %   as regular layers.
    %
    %   TO 110810; 120516
    %
    %   Now, eind 2015 with years of experience, I decided that confining
    %   beds have zero thickness, as in old-fashioned quasi-3D models that
    %   had both zero thickness aquifers (with tranmissivity) and zero
    %   thickness aquitards (with resistance or its inverse).
    %   For years confing beds as in 3D models has seldom been useful and
    %   pratically always confing beds were defined as regular model
    %   layers. An old fashioned confing bed may be useful to simulate
    %   for instance a low-conductivity holocene coverlayer as it allows
    %   to concentrate its resistance at its bottom. Notice, however that
    %   you can achieve the same through the vertical conductance of the cover
    %   layer.
    %
    %   TO 151115
    
    %% Constant properties
    properties (Constant) % also physically stored
    end
    %% Directly accessible properties
    properties
        %%
        % Properties that are actually stored. Their direct setting is
        % discouraged as this may cause inconsistencies. Always use the
        % standard constructor call
        % gr = gridObj(xGr,yGr,zGr,'LAYCBD',LAYCBD,'MINDZ',MINDZ,'AXIAL',true|false);
        %
        MINDZ=0.001;             % minimum layer thickness
        AXIAL=false;             % ~=0 if axial symmetric flow
        LAYCBD=0;                % layer confining bed vector
        isLay;                   % vector Nlay+Ncbd long telling if layer is LAY or CBD
        xGr; yGr,zGr             % grid line coordinates, vectors of size Nx+1, Ny+1, Nz+1
        eGr,nGr                  % google earth grid coordinates, easting and northing degrees
        latlon                   % boolean true of lonGr and latGr are given
        xw0=0; yw0=0; zw0=0; anglew=0;  % position of model zero in real world coordinates
        layersAreUniform,        % true of all layers are uniform
        UserData;                % for any purpose
    end
    %% Protected properties
    properties (Access=protected)
        %%
        % To guarantee integrety of the grid this grid is protected for
        % direct access by the user. It can only be set by the constructor
        % and by methods of the object
        Zfull
    end
    %% Dependent properties
    properties (Dependent=true)
        %%
        % Properties that are computed as and when needed.
        size, sizeLay, sizeCBD, sizeFULL, % 3rd dimension size resepectively: Nlay Nlay Ncbd Nz
        xlim,ylim,
        zlim,  % applied to zGr,in fact only useful for uniform grid
        
        Nxy, Nxyz   % Nx*Ny,   Nx*Ny*Nlay (total nr of model cells)
        
        %%
        % Indices of layers and confining bed in total layer stack
        ITlay,ITcbd,IBlay,IBcbd, % pointer to top of layer  LAYTOP=Z(:,:,Itop)
        
        %%
        % Size of grid
        Nx, Ncol, Ny, Nrow, Nz,Nlay, Ncbd,    % number of cells in direction
        
        %%
        % Coordinates
        % xm,xh=xc are cell centers. xh,xc same as xm with outer coordinates
        %  replaced by xGr([1 end]) for contouring. xp is xGr(2:end-1) for
        %  contouring stream function. Xm is full 2D cell center
        %  coordinates. Same for y and z
        % Variables that start with a capital letter are 3D arrays (except
        % Xm and Ym which are 2D arrays) while varibles starting with a
        % lower case letter are vectors.
        xm, Xm, XM, XMlay, xh, xc,Xc, XC, xp, XP, Xp, XGR, XGr,
        ym, Ym, YM, YMlay, yh, yc,Yc, YC, yp, YP, Yp, YGR, YGr,          
        xkm, ykm   % distance in km
%         em, ec   % same as xm xc but in wgs longitude degrees
%         nm, nc   % same as ym yx but in wgs latitude  degrees                
        
        ZGR, Z, zh, zc, ZC, zp, ZP, Zp, ZMlay, ZMcbd, ZM,
        zlay, zLay, zTlay, zmlay, zMlay, zBlay,
        zcbd,       zTcbd, zmcbd, zMcbd, zBcbd,
        ZGrTlay, ZGrBlay, ZGrTcbd, ZGrBcbd,  % full grid of tops and bottoms of Lay and CBD
        
        Zlay, ZGRlay, ZTlay, ZBlay, ZClay, ZTgrlay, ZMgrlay, ZBgrlay,
        Zcbd, ZGRcbd, ZTcbd, ZBcbd, ZCcbd,
        
        %%
        % Distances and thicnesses
        dx,           DXlay, DXcbd, DX,
        dy,           DYlay, DYcbd, DY,
        dzlay, dzcbd, DZlay, DZcbd, DZ,
        dz, zm,
        
        %%
        % Volumes and aeas
        Vlay, vlay             % volume of model cells and of entire model
        Vcbd, vcbd             % volume of CBD per cell and total
        AREA, area             % surface area of model cells and entire model
        AREA3,                 % surface area for all cells in layers (Ny,Nx,Nlay)
        %%
        % distance to xm=0 (along the x-axis) for axi-symmetric situations
        r, rm, dr, R, RM, rkm
        
        %%
        % Axial flow
        TWOPIR,     % 2*pi*(abs(XM)).^2, full size
end
    methods
        %% Constructor
        function o=gridObj(xGr,yGr,Z,varargin)
            %% gr = gridObj(xGr,yGr,Z,'LAYCBD',LAYCBD,'MINDZ',MINDZ,'AXIAL',true|false)
            % Constructor for the gridObj.
            % Note that there is no way to know that the Z-array associates
            % correctly with the xGr and the yGr. This is the user's
            % responsibility.
            %
            % TO 151115
            
            if nargin==0;
                o=gridObj([0 1],[1 0],[1 0]);
                return;
            end
            
            if nargin<3,
                error('not enough input arguments for gridObj: use gridObj(xGr,yGr,zGR[,LAYCBD[,MINDX[,AXIAL]]])');
            end

            if isempty(xGr), error('%s: Empty first argument xGr not allowed.',mfilename); end
            if isempty( Z ), error('%s: Empty third argument zGr|Z not allowed.',mfilename); end            
            if isempty(yGr)  % may be used for flat models (1 row wide)
                yGr = [ -0.5 0.5];
            end

            % allow for variablNm,value pair input like
            % gridObj(xGr,yGr,zGr,'AXIAL',true);
            [o.latlon,varargin] = getProp(varargin,{'latlon','ll','lat','wgs'},o.latlon);
            [o.LAYCBD,varargin] = getProp(varargin,'LAYCBD',o.LAYCBD);
            [o.MINDZ ,varargin] = getProp(varargin,'MINDZ' ,o.MINDZ);
            [o.AXIAL ,varargin] = getProp(varargin,'AXIAL' ,o.AXIAL);
            
            % but also regular input like  gridObj(xGr,yGr,zGr,LAYCBD,MINDZ,AXIAL)
            [o.LAYCBD,varargin] = getNext(varargin,{'logical','double'},o.LAYCBD);
            [o.MINDZ ,varargin] = getNext(varargin,'double',o.MINDZ);
            [o.AXIAL ,varargin] = getNext(varargin,{'logical','double'},o.AXIAL);
            
            if ~(islogical(o.AXIAL) || (isnumeric(o.AXIAL) && isscalar(o.AXIAL)))
                error('%s: AXIAL must be a logical or numeric scalar',mfilename);
            end
            
            if ~isempty(varargin)
                msgId = sprintf('mfLab:%s:vararginNotUsed',mfilename);
                warning('on',msgId);
                warning(msgId,'%s: varargin not fullY used',mfilename);
            end
            
            % We strictly adhere to the node numbering in MODFLOW. That is,
            % the x-axis is increasing from west to east.
            % The y axis is decreasing from north to south, while the
            % north-west corner of the model has indices 1.1.
            % The z axis is decreasing from top to bottom. The top layer
            % being layer number one.
            % NW corner of the model
            
            % if latlon option is used, xGr and yGr area assumed to be in
            % lon, lat respectively.
            if o.latlon
                o.eGr = xGr;
                o.nGr = yGr;
                [xGr,  ~] = wgs2utm(o.eGr,mean(o.nGr)*ones(size(xGr)));
                [~  ,yGr] = wgs2utm(mean(o.eGr)*ones(size(yGr)),o.nGr);
            end
            
            o.xGr = unique(roundn(xGr(:),4))'; % Guarantee x is unique and increasing
            
            o.yGr = unique(roundn(yGr(:),4));
            %%
            % Enforce y to be decreasing from north to south
            o.yGr=o.yGr(end:-1:1); % Guarantee y is decreasing
            
            %%
            % Z can be a full 3D array with UL corner corresponding to NW
            % top of model, where Z has the top and the bottom of all
            % layers (LAY and CBD) in sequence for each cells, i.e. not the
            % cell corners. Cell corners can be computed as the averate of
            % the cell elevations of each pair of adjacent cells. This is
            % done in the ZGR method of this gridObj.
            % Z can also be a vector with the layer elevations. This is the
            % case of all layers are uniform. In that case the gridObj only
            % stores this vector, while the methods of this gridObj can
            % still compute full 3D versions if requested.
                        
            if numel(size(Z))==3 && size(Z,1)==1 && size(Z,2)==1
                % Z is a vector in the third dimension
                Z=XS(unique(roundn(XS(Z),4))); % treat Z as a vector like zGr
                % force Z to decrease
                o.layersAreUniform=true;
                o.Zfull = Z(end:-1:1);
                o.zGr   = o.Zfull;
            elseif isvector(Z)
                o.layersAreUniform=true;
                Z = XS(unique(roundn(XS(Z(:)),4)));
                o.Zfull = Z(end:-1:1);
                o.zGr   = o.Zfull;
            elseif size(Z,3)==1, % Z must be 3D, i.e. Z must have at least 2 values in Z-direction.
                    error('gridObj:gridObj:Znot3D',...
                        '3rd argument Z (or zGr) must be 3D unless it is given as a vector');
            else
                %%
                % At this point we have a full 3D Z array.
                % We define for this situation o.zGr as the vector
                % of layer elevation averages.
                % Even in this case having this the vector o.zGr may be usefull
                % in some cases, although it is not compatible with a full-size Z
                % in which each column may have different elevations.
                
                %%
                % first check the x,y dimensions of Z
                if size(Z,1) ~= size(o.yGr,1)-1
                    error('%s: y-dimension of Z (%d) does not match length(yGr)-1 (%d)',...
                        mfilename,size(Z,1),size(o.yGr,1)-1);
                end
                if size(Z,2) ~= size(o.xGr,2)-1
                    error('%s: x-dimension of Z (%d) does not match length(xGr)-1 (%d)',...
                        mfilename,size(Z,2),size(o.xGr,2)-1);
                end
                
                o.zGr = roundn(mean(mean(Z,1),2),4);
                
                %%
                % Even if we have a full 3D Z-array layers may all be
                % uniform. So check it
                o.layersAreUniform = all(Z(1,1,:)==o.zGr);
                if o.layersAreUniform
                    %%
                    % If layers are all uniform, then o.zGr contains all
                    % non-redundant grid information and this is all we
                    % need to store.
                    o.zGr   = XS(unique(YS(o.zGr)));
                    o.zGr   = o.zGr(end:-1:1);
                    o.Zfull = o.zGr;
                else % we have to store the full 3D Z-array
                    o.Zfull = roundn(Z,4);
                    if diff(o.zGr([1 end]))>0
                        o.zGr   = flip(o.zGr,3);
                        o.Zfull = flip( Z    ,3);
                    else
                        o.Zfull = Z;
                    end
                end
            end

            %% guarantee size(..,i) is >=2 in all three i directions
            % The grid coordinate vectors and or arrays must have at
            % least two values in every direction
            %%
            % x-directon
            if numel(o.xGr)==1 || ~isvector(o.xGr),
                error('gridObj/gridObj: arg 1, xGr must be a vector with >1 elements');
            end
            % y-direction
            if numel(o.yGr)==1 || ~isvector(o.yGr),
                error('gridObj/gridObj: arg 2, yGr, must be a vector with >1 elmements.');
            end
            % z-direction (has already been verified above.
            
            %% Make sure Z runs from high to low
            % We already made sure above that x is increasing and y is
            % decreasing. But we did not yet check that z is decreasing.
            %
            % Perhaps Z does not exist over some parts of the model area
            % due to inactive cells. To check the direction of the grid
            % (increasing or decreasing) we must have at least one z-column
            % with non-NaN values at top and bottom. This is not trivial at
            % all in a general meesh, especially if layer thicknesses are
            % allowed to be zero.            
            %%
            % Find cells in top layer with non-NaN elevation
            I=find(~isnan(o.Z(:,:,1)));
            %%
            % Find the first corresponding non-NaN cell cell of the bottom layer
            for ib=1:length(I)
                
                % Global index of corresponding bottom cell
                idx = I(ib)+o.Ny*o.Nx*(size(o.Z,3)-1);
                
                if ~isnan(o.Z(idx))
                    break;
                elseif ib == length(I)
                    error('gridObj:gidObj:zBotNaN',...
                        'There are no combinations of top and bottom cells where both are non-NaN');
                end
            end
            
            %%
            % Issue an error if the grid is non vertically decreasing. We
            % are not going to solve it as it may violate the integrety of
            % the grid without the user knowing it.
            if o.Z(I(ib))<o.Z(idx)
                error('gridObj:gridObj:ZnotVerticallyDecreasing',...
                    '%s Z (elevation) array must be decreasing along the 3rd dimension',mfilename);
            end

            %% Guarantee minimum layer thickness            
            minLayerThickness=roundn(min(min(min(-diff(o.Zfull,1,3)))),3);
            if minLayerThickness<roundn(o.MINDZ,3)
                msgId = 'gridOb:DZsmallerMinDZ';
                
                warning('on',msgId);
                warning(msgId,['Minimum layer thickness =%g < MINDZ = %g\n'...
                    'Generating spy figure to show in which layers and where this happens.'],...
                    minLayerThickness,o.MINDZ);
                warning('off',msgId);
                
                o.checkDZ;
                error(['%s: Minimum layer thickness = %g (<MINDZ = %g).\n',...
                     'Remedy: Correct this before restarting.\n'],...
                     mfilename,minLayerThickness,o.MINDZ);                 
            end

            %% Guarantee Z array of size Ny,Nx,Nz+1
            % If Z is of size Ny+1,Nx+1,Nz+1 make if of size Ny,Nz,Nz+1
            if ~o.layersAreUniform
                if size(o.Zfull,1)>o.Ny,  o.Zfull=0.5*(o.Zfull(1:end-1,:,:)+o.Zfull(2:end,:,:)); end
                if size(o.Zfull,2)>o.Nx,  o.Zfull=0.5*(o.Zfull(:,1:end-1,:)+o.Zfull(:,2:end,:)); end
                o.zGr = mean(mean(o.Zfull,1),2);
            end
            
            %% The grid may represent an axially symmetric model
            % Axial symmetry is always along the x-axis then representing
            % the r-axis. We want the rGr-axis (xGr-axis) to containt coordinate
            % zero. If necessary insert xGr=0 in both xGr and zGr arrays
            if o.AXIAL
                i=find(o.xGr<0,1,'last');
                j=find(o.xGr>0,1,'first');
                if ~isempty(i) && ~isempty(j) && isempty(o.xGr==0)
                    
                    % Insert xGr=0
                    o.xGr=[o.xGr(1:i  )         1            o.xGr(j:end)];
                    % also split Z (but only if Z is a full 3D array)
                    if ~o.layersAreUniform
                        o.Zfull    =[o.Zfull(:,1:i,:) mean(o.Zfull(:,i:j,:),2) o.Zfull(:,j:end,:)];
                    end
                end
            end
            
            o.LAYCBD=isLayer(o.Nz,o.LAYCBD);
                                    
        end
        
        function Id = IdTop(o,Idx)
            %GRIDOBJ/IDTOP computes the indices within the top layer of the
            %grid given Idx for any layer in the grid
            % TO 131011?
            Id = rem(Idx-1,o.Nxy)+1;
        end
        
        function o = convert(o,Ix,Iy)
            %% Convert the grid to reflect then new outcut defined by Ix,Iy
            %  where Ix and Iy are the old indices.
            if max(Ix)>o.Nx || min(Ix)<1,
                error('gridObj/convert: Ix indices must be between 1 and %d.',o.Nx); end
            if max(Iy)>o.Ny || min(Iy)<1,
                error('gridObj/convert: Iy indices must be between 1 and %d.',o.Ny); end
            
            o.Zfull  = o.Zfull(Iy,Ix,:);
            o.xGr    = o.xGr([ Ix(1:end), Ix(end)+1 ]);
            o.yGr    = o.yGr([ Iy(1:end); Iy(end)+1 ]);

        end
        
        function AXIAL = get.AXIAL(o), AXIAL = o.AXIAL;  end
        function LAYCBD= get.LAYCBD(o)
            LAYCBD=o.LAYCBD; end
        
        function size    = get.size(o),     size     = [o.Ny, o.Nx, o.Nlay]; end
        function sizeLay = get.sizeLay(o),  sizeLay  = [o.Ny, o.Nx, o.Nlay]; end
        function sizeCBD = get.sizeCBD(o),  sizeCBD  = [o.Ny, o.Nx, o.Ncbd]; end
        function sizeFULL= get.sizeFULL(o), sizeFULL = [o.Ny, o.Nx, o.Nz  ]; end
        
        function Nxy  = get.Nxy(o),  Nxy = o.Nx*o.Ny;         end
        function Nxyz = get.Nxyz(o), Nxyz= o.Nx*o.Ny*o.Nlay;  end
        function Nx   = get.Nx(o) ,  Nx   = numel(o.xGr)-1;   end
        function Ncol = get.Ncol(o), Ncol = numel(o.xGr)-1;   end
        function Ny   = get.Ny(o),   Ny   = numel(o.yGr)-1;   end
        function Nrow = get.Nrow(o), Nrow = numel(o.yGr)-1;   end
        function Nz   = get.Nz(o) ,  Nz   = numel(o.zGr)-1;   end
        function Nlay = get.Nlay(o), Nlay = o.Nz;             end
        function Ncbd = get.Ncbd(o), Ncbd = sum(o.LAYCBD > 0);end
        
        %% X coordinates
        function x    = xBox(o), x = o.xGr([1 end end 1   1]); end
        function y    = yBox(o), y = o.yGr([1 1   end end 1]); end
        function xlim=  get.xlim(o), xlim = [min(o.xGr) max(o.xGr)]; end
        function xGr  = get.xGr(o), xGr = o.xGr; end
        function xm   = get.xm(o),  xm  = 0.5*(o.xGr(1:end-1)+o.xGr(2:end));  end
        function xkm  = get.xkm(o), xkm  = o.xm/1000; end        
        function XGR  = get.XGR(o), XGR = bsxfun(@times,o.xGr,ones(o.Ny+1,1,o.Nz+1)); end
        function XGr  = get.XGr(o), XGr = bsxfun(@times,o.xGr,ones(o.Ny+1,1,1)); end
        function XM   = get.XM(o),  XM  = bsxfun(@times,o.xm ,ones(o.Ny,1,o.Nz)); end
        function XMlay= get.XMlay(o),XMlay = o.XM; end        
        function Xm   = get.Xm(o),  Xm  = bsxfun(@times,o.xm ,ones(o.Ny,1,1)); end        
        function Xc   = get.Xc(o)
            %GRIDOBJ/YC -- geneates coordinates for plotting in xz plane,
            %outer coordinates are on xGr inner ones on xm
            %
            % USAGE: Xc = gr.Xc();
            %    size of Xc if [gr.Ny,gr.Nx]

            Xc = o.Xm;
            if o.Nx>1,
                Xc(:,  1)=o.xGr(  1);
                Xc(:,end)=o.xGr(end);
            end
        end
        function XC   = get.XC(o)
            %GRIDOBJ/XC -- geneates coordinates for plotting in xz plane,
            %outer coordinates are on xGr inner ones on xm
            %
            % USAGE: Xc = gr.Xc();
            %    size of Xc if [gr.Ny,gr.Nx,gr.Nz]
            XC = o.XM;
            if o.Nx>1,
                XC(:,  1,:)=o.xGr(  1);
                XC(:,end,:)=o.xGr(end);
            end
        end

        function XP = get.XP(o)
            %GRIDOBJ/XP -- geneates x coordinates for plotting the stream function (xz-planes)
            %
            % USAGE: XP = gr.XP()
            % size of XP is [gr.Nx-1,gr.Nz+1]            
            XP  = bsxfun(@times,o.xGr(1,2:end-1),ones(1,1,o.Nz+1)); % always use 1st xz plane (first iy==1);
        end


        function ylim=  get.ylim(o), ylim = [min(o.yGr) max(o.yGr)]; end
        function yGr   = get.yGr(o), yGr = o.yGr; end
        function ym    = get.ym( o), ym  = 0.5*(o.yGr(1:end-1)+o.yGr(2:end));  end
        function ykm   = get.ykm(o), ykm  = o.ym/1000; end        

        function YGR   = get.YGR(o),  YGR = bsxfun(@times,o.yGr,ones(1,o.Nx+1,o.Nz+1)); end
        function YGr   = get.YGr(o),  YGr = bsxfun(@times,o.yGr,ones(1,o.Nx+1,1)); end
        function YM    = get.YM(o),   YM  = bsxfun(@times,o.ym,ones(1,o.Nx,o.Nz)); end
        function YMlay = get.YMlay(o),YMlay = o.YM; end        
        function Ym    = get.Ym(o),   Ym  = bsxfun(@times,o.ym,ones(1,o.Nx,1));    end
        function YC    = get.YC(o)
            %GRIDOBJ/YC -- geneates coordinates for plotting in xz plane,
            %outer coordinates are on yGr inner ones on ym
            %
            % USAGE: YC = gr.YC();
            %    size of YC if [gr.Ny,gr.Nx,gr.Nz]
            YC = o.YM;
            if o.Ny>1
                YC(1,  :,:)=o.yGr(  1);
                YC(end,:,:)=o.yGr(end);
            end
        end
        function Yc    = get.Yc(o)
            %GRIDOBJ/YC -- geneates coordinates for plotting in xy plane,
            %outer coordinates are on yGr inner ones on ym
            %
            % USAGE: Yc = gr.Yc();
            %    size of Yc if [gr.Ny,gr.Nx]
            
            Yc = o.Ym;
            if o.Ny>1
                Yc(1,  :)=o.yGr(  1);
                Yc(end,:)=o.yGr(end);
            end
        end
        function YP = get.YP(o)
            %GRIDOBJ/YP -- generates coordinates for plotting stream function (xz-planes)
            %
            % USAGE: YP = gr.YP();
            % size of YP is [gr.Ny-1,gr.Nz+1]            
            YP  = bsxfun(@times,o.yGr(2:end-1),ones(1,1,o.Nz+1)); % always use 1st xz plane (first iy==1);
        end


        %% Z-coordinates
        function zlim  = get.zlim(o)  , zlim  = [min(o.zGr) max(o.zGr)]; end
        function zm    = get.zm(o)    , zm    = 0.5*(o.zGr(1:end-1)+o.zGr(2:end));  end  % center of model layers
        
        function zmlay = get.zmlay(o) , zmlay = o.zm;  end  % center of model layers
        function zmcbd = get.zmcbd(o) , zmcbd = o.zm(o.LAYCBD>0);  end  % center of model layers        

        function zMlay = get.zMlay(o) , zMlay = o.zmlay;  end  % center of model layers
        function zMcbd = get.zMcbd(o) , zMcbd = o.zmcbd;  end  % center of model layers

        function zlay  = get.zlay(o),zlay   = o.zGr; end
        function zLay  = get.zLay(o),zLay   = o.zlay; end
        function zTlay = get.zTlay(o),zTlay = o.zGr(2:end);   end
        function zBlay = get.zBlay(o),zBlay = o.zGr(2:end); end
        function zTcbd = get.zTcbd(o),zTcbd = o.zGr(     o.LAYCBD>0   ); end
        function zBcbd = get.zBcbd(o),zBcbd = o.zGr(find(o.LAYCBD>0)+1); end
                
        function ZGR = get.ZGR(o)
            %% All nodes of the mesh !!
            if o.layersAreUniform
                ZGR = bsxfun(@times,o.zGr,ones(o.Ny+1,o.Nx+1,1));
            else
                ZGR = 0.5*(cat(2,o.Zfull(:,1,:),o.Zfull) + cat(2,o.Zfull,o.Zfull(:,end,:)));
                ZGR = 0.5*(  ZGR  ([1 1:end],:,:) +   ZGR  ([1:end end],:,:));
            end
        end
                                
        function Z   = get.Z(o)
            if o.layersAreUniform, Z = bsxfun(@times,ones(o.Ny,o.Nx,1),o.zGr);
            else                   Z = o.Zfull;          end 
        end
        function ZM    = get.ZM(o), ZM  = 0.5*(o.Z(:,:,1:end-1)+o.Z(:,:,2:end)); end
        function ZC    = get.ZC(o) 
            ZC  = o.ZM;
            if o.Nz>1
                ZC(:,:,[1 end])= o.Z(:,:,[1 end]);
            end
        end
        function ZP    = get.ZP(o)  % for plotting stream function (xz-planes)
            if size(o.Z,2)>1
                ZP  = 0.5*(o.Z(:,1:end-1,:)+o.Z(:,2:end,:));
            end
            if size(ZP,1)>1
                ZP  = 0.5*(o.Z(1:end-1,:,:)+o.Z(2:end,:,:));
            end
        end
        function Zlay   = get.Zlay( o), Zlay = o.Z; end
        function ZTlay  = get.ZTlay(o), ZTlay= o.Z(:,:,1:end-1); end
        %%
        % Z-values at xGr,yGr
        function ZTgrlay = get.ZTgrlay(o), ZTgrlay = o.ZGR(:,:,1:end-1); end
        function ZBgrlay = get.ZBgrlay(o), ZBgrlay = o.ZGR(:,:,2:end  ); end
        function ZMgrlay = get.ZMgrlay(o), ZMgrlay = 0.5*(o.ZTgrlay + o.ZBgrlay); end
        
        function ZMlay = get.ZMlay(o), ZMlay = o.ZM; end  % center Lay+CBD
        function ZBlay = get.ZBlay(o), ZBlay = o.Z(:,:,2:end); end
        function ZClay = get.ZClay(o), ZClay = o.ZM;
            if o.Nlay>1
                ZClay(:,:,  1) = o.ZTlay(:,:,  1);
                ZClay(:,:,end) = o.ZBlay(:,:,end);
            end
        end

        function Zcbd= get.Zcbd(o),Zcbd = o.Z( :,:,o.LAYCBD>0  ); end
        function ZTcbd = get.ZTcbd(o), ZTcbd = o.Zcbd; end 
        function ZMcbd = get.ZMcbd(o), ZMcbd = o.Zcbd; end
        function ZBcbd = get.ZBcbd(o), ZBcbd = o.Zcbd; end 
   
        function dx     = get.dx(o),    dx     =  abs(diff(o.xGr,1,2)); end
        function dy     = get.dy(o),    dy     =  abs(diff(o.yGr,1,1)); end
        function dz     = get.dz(o),    dz     =  abs(diff(o.zGr,1,3)); end
        function dzlay  = get.dzlay(o), dzlay  = o.dz;    end
        function dzcbd  = get.dzcbd(o), dzcbd  = o.dz(:,:,o.LAYCBD>0) * 0; end

        function DX     = get.DX(o),  DX     = bsxfun(@times,ones(o.Ny,1,o.Nz),o.dx); end
        function DXlay  = get.DXlay(o), DXlay  = o.DX; end
        function DXcbd  = get.DXcbd(o), DXcbd  = o.DX(:,:,o.LAYCBD>0); end
        function DY     = get.DY(o),  DY     = bsxfun(@times,ones(1,o.Nx,o.Nz),o.dy); end
        function DYlay  = get.DYlay(o), DYlay  = o.DY; end
        function DYcbd  = get.DYcbd(o), DYcbd  = o.DY(:,:,o.LAYCBD>0); end
        function DZ     = get.DZ(o)
            if isvector(o.Z)
                DZ = bsxfun(@times,ones(o.Ny,o.Nx,1),-diff(o.Z,1,3));
            else
                DZ = abs(diff(-o.Z,1,3));
            end
        end
        function DZlay  = get.DZlay(o), DZlay  = o.DZ;  end
        function DZcbd  = get.DZcbd(o), DZcbd  = o.DZ(:,:,o.LAYCBD>0) * 0; end
        
        function Vlay = get.Vlay(o),...
            if o.AXIAL
                Vlay = pi*(o.XGR(1:end-1,2:end,:).^2-o.XGR(1:end-1,1:end-1,:).^2).* o.DZlay;
            else
                Vlay = o.DXlay.* o.DYlay.* o.DZlay;
            end
        end
        function vlay = get.vlay(o), vlay = sum(o.Vlay(:));           end
        function Vcbd = get.Vcbd(o)
            if o.AXIAL
                Vcbd = pi*(o.XGR(1:end-1,2:end,o.LAYCBD>0).^2-o.XGR(1:end-1,1:end-1,o.LAYCBD>0).^2).* o.DZcbd;
            else
                Vcbd = o.DXcbd.* o.DYcbd.* o.DZcbd;
            end
        end
        function vcbd = get.vcbd(o); vcbd = sum(o.Vcbd(:));           end
        
        function AREA3= get.AREA3(o), AREA3= bsxfun(@times,o.AREA,ones(1,1,o.Nlay)); end
        function AREA = get.AREA(o)
            if o.AXIAL
                AREA = ones(size(o.dy)) * (pi.*(o.xGr(2:end).^2-o.xGr(1:end-1).^2));
            else
                AREA = o.dy * o.dx;    
            end
        end
        function area = get.area(o),  area = sum(sum(o.AREA)); end

        % axial symmetric, assuming r along x-axis using abs(x)
        function r  = get.r(o)
            % GRIDOBJ/R -- gets radius from x=0
            %
            % USAGE: gr.r()
            
            if o.AXIAL
                if o.Nx==1
                    r  = sqrt(o.yGr.^2);
                else
                    r  = sqrt(o.xGr.^2);
                end
            else
                r = sqrt(bsxfun(@times,o.yGr,ones(1,o.Nx+1)).^2 + ...
                         bsxfun(@times,ones(o.Ny+1,1),o.xGr).^2);
            end
        end
        function rm = get.rm(o)
            if o.AXIAL
                if o.Nx==1
                    rm = abs(o.ym);
                else
                    rm = abs(o.xm);
                end
            else
                rm = sqrt(bsxfun(@times,o.ym,ones(1,o.Nx)).^2 + ...
                         bsxfun(@times,ones(o.Ny,1),o.xm).^2);
            end
        end
        function rkm = get.rkm(o), rkm = o.rm/1000; end
        function dr = get.dr(o)
            if o.AXIAL
                if o.Nx==1
                    dr = abs(o.dy);
                else
                    dr = abs(o.dx);
                end
            else
                dr = NaN(o.Ny,o.Nx);
            end
            
        end
        function R  = get.R( o)
            if o.AXIAL
                if o.Nx==1
                    R  = sqrt(o.YGR.^2); 
                else
                    R  = sqrt(o.XGR.^2); 
                end
            else
                R = sqrt(o.XGR.^2+o.YGR.^2);
            end
        end
        function RM = get.RM(o)
            if o.AXIAL
                if o.Nx==1
                    RM = sqrt(o.YM.^2);
                    % tackle cases where center node is at x==0, prevent
                    % that we have r==0 at this locattion, use 
                    I=find(RM==0);
                    if ~isempty(I)
                        RM(I) = sqrt(0.5*o.YM(I+1).^2); % xGr to the right of xm(I) => xGr(I+2)
                    end
                else
                    RM = sqrt(o.XM.^2);
                    % tackle cases where center node is at x==0, prevent
                    % that we have r==0 at this locattion, use 
                    I=find(RM==0);
                    if ~isempty(I)
                        RM(I) = sqrt(0.5*o.XM(I+1).^2); % xGr to the right of xm(I) => xGr(I+2)
                    end
                end
            else
                RM = sqrt(o.XM.^2 + o.YM.^2);
            end
        end
 
        function TWOPIR = get.TWOPIR(o), TWOPIR=2*pi*o.RM; end
        
        % world coordinates
        function xw0     = get.xw0(o), xw0=o.xw0; end
        function yw0     = get.yw0(o), yw0=o.yw0; end
        function zw0     = get.zw0(o), zw0=o.zw0; end
        function anglew  = get.anglew(o), anglew=o.anglew; end
        
        % facilitates plotting xh,yh,zh phased out use xc,yc,hc
        function xh  = get.xh(o) , xh  =  o.xc; end  % for heads
        function yh  = get.yh(o) , yh  =  o.yc; end  % for heads in xy plane
        function zh  = get.zh(o) , zh  =  o.zc; end  % for heads in zx plane

        function xc  = get.xc(o) , xc  = o.xm;  xc([1 end])=   o.xGr([1 end]);     end  % for concentrations  
        function yc  = get.yc(o) , yc  = o.ym;  yc([1 end])=   o.yGr([1 end]);     end  % for concentrations  
        function zc  = get.zc(o) , zc  = XS(o.zm); zc([1 end])=XS(o.zGr([1 end])); end  % for concentrations
        
        % coordinates to plot Psi (on all horizontal grid lines and
        % vertical gridline except the outer two
        function xp  = get.xp(o) , xp= o.xGr(2:end-1); end  % for Psi
        function yp  = get.yp(o) , yp= o.yGr(2:end-1); end  % for Psi
        function zp  = get.zp(o) , zp= XS(o.zGr);      end  % for Psi
        function Xp  = get.Xp(o) , Xp= o.XP; end
        function Yp  = get.Yp(o) , Yp= o.YP; end
        function Zp  = get.Zp(o) , Zp= o.ZP; end
        
        % index in array Z which combines all interfaces
        function ITlay = get.ITlay(o), ITlay = 1:o.Nz  ; end
        function IBlay = get.IBlay(o), IBlay = 2:o.Nz+1; end
        function ITcbd = get.ITcbd(o), ITcbd = find(o.LAYCBD>0)  ; end
        function IBcbd = get.IBcbd(o), IBcbd = find(o.LAYCBD>0)+1; end

        function const = const(o,value)
            %GRIDOBJ/CONST -- generates a const array of size gr.size
            %
            % USAGE: array = gridObj.const(value|vector)
            %
            % If value is scalar:
            %   generates a 3D array of size [gr.Ny,gr.Nx,gr.Nlay], with all values=value.
            % If value is vector:
            %   generates a 3D array of size [gr.Ny,gr.Nx,length(value)], with
            %   array(:,:,i)=value(i);
            %
            % TO 120410
            if isempty(value) || ~(isnumeric(value) || islogical(value))
                error('%s: input must be numeric and not empty',mfilename);
            elseif numel(size(value)) == numel(o.size)  &&  all(size(value)==o.size)
                % do nothing
                const = value;
            elseif numel(size(value))==2 && all(size(value)==[o.Ny,o.Nx])
                % expand plane in z-direction
                const = bsxfun(@times,ones(1,1,o.Nlay),value);
            elseif isscalar(value)
                % expand scalar in z-direction
                const=ones(o.Ny,o.Nx,o.Nlay)*value;
            else % assume vector is a z-vector
                % expand over plance. Keep length of z in z-direction
                const=bsxfun(@times,ones(o.Ny,o.Nx,1),XS(value(:)));
            end
            if islogical(value)
                const = logical(const);
            end
        end
        
        function const = constLay(o,value)
            %GRIDOBJ/CONSTLAY -- generate const array of size of layers of grid
            %
            % USAGE: const = gr.constLay(value);
            %
            %LAY same as const = gr.const(value)
            % TO 120101
            
            const = o.const(value);
        end
        
        function const = constCB(o,value)
            %GRIDOBJ/CONSTCB -- generates an array for all confining beds (CB)
            %
            % USAGE: array = gridObj/constCB(value or vector)
            %
            % If value is scalar:
            %    generates a 3D array of size [gr.Ny,gr.Nx,gr.Ncbd], with all values=value.
            % If value is vector:
            %   generates a 3D array of size [gr.Ny,gr.Nx,length(value)], with
            %   array(:,:,i)=value(i);
            %
            % TO 130217
            if isscalar(value)
                const=ones(o.Ny,o.Nx,o.Ncbd)*value;
            else
%                 if numel(value)~=o.Ncbd
%                     error('%s/array(..) requires one argument, a scalar or a vector length Nz\n',mfilename);
%                 end
                const=bsxfun(@times,ones(o.Ny,o.Nx,1),XS(value(:)));
            end
        end
        
        function ZGrTlay = get.ZGrTlay(o) % ZGR of layer tops
            ZGrTlay = o.ZGR(:,:,1:o.Nz);
        end
        function ZGrBlay = get.ZGrBlay(o) % ZGR of layer bottoms
            ZGrBlay = o.ZGR(:,:,2:o.Nz+1);
        end
        function ZGrTcbd = get.ZGrTcbd(o)% ZGR of cbd tops
            ZGrTcbd = o.ZGR(:,:,o.LAYCBD>0);
        end
        function ZGrBcbd = get.ZGrBcbd(o)% ZGR of cbd tops
            ZGrBcbd = o.ZGR(:,:,(o.LAYCBD>0) +1);
        end
        function checkDZ(o)
            %GRIDOBJ/CHECKDZ -- shows using spy where thickness of layer is < MINDZ
            %
            % USAGE gr.checkDZ()
            %
            % TO 130925
            
            figure('name',sprintf('Checking if layer thickness < %g',o.MINDZ),'pos',screenPos(0.75));
            
            n = ceil(sqrt(size(o.DZ,3)));
            for i=1:size(o.DZ,3)
                subplot(n,n,i);
                spy(o.DZ(:,:,i)<o.MINDZ);
                title(sprintf('DZ(%d)<%g?',i,o.MINDZ));
            end
        end
        function dist = dist(o,x,y)
            %GRIDOBJ/DIST -- yields the distance of all cell centers to point x,y
            %
            % USAGE: gridObj/dist: dist=gr.dist(x,y)
            %
            % TO 120410
            
            switch nargin
                 case 1,
                     x=0; y=0;
                 case 3,
                     if ~isnumeric(x) || ~isnumeric(y)
                         error('gridObj/dist: x and y must be scarlars.');
                     end
                 otherwise
                         error('gridObj/dist: wrong number of argumens gridObj.dist() or gridObj.dis(x,y)');
            end
            dist=sqrt((o.XMlay-x).^2+(o.YMlay-y).^2);
        end
        function o=setWorld(o,xw0,yw0,zw0,angle)
            %GRIDOBJ/SETWORLD -- sets world coordinate system for this grid
            %
            % USAGE: gridObj/setWorld: obj=gr.setWorld(xw0,yw0,anglew)
            %
            % The coordinates xw0,yw0 match model coordinates 0,0
            % angle is rotation of model relative to world, anticlockwise
            % from hoirzontal (east) as in mathematics.
            %
            % TO 120420
            
            if nargin<3, error('gridObj:setWorld:notEnoughInputArgs',...
                    ['gridObj/setWorld: insuffcient input arguments, use\n',...
                     '(xw,yw) or (xw,yw,anglew) or (xw,yw,zw,anglew).']);
            else
                o.xw0=xw0;
                o.yw0=yw0;
                if nargin==4
                    o.anglew=zw0;
                else
                    o.zw0=zw0;
                    o.anglew=angle;
                end
            end
            o.anglew = o.anglew*pi/180; % relative to east (x-axis)
        end
        
        function [xw, yw]=world(o,xm,ym)
            %GRIDOBJ/WORLD -- computes xw yw in world coordinates
            %
            % USAGE: gridObj/world:  [xw,yw]=gr.world(xm,ym)
            %
            % computes xw yw in world coordinates, from model coordinates xm,ym
            % make sure you use consistent coordiates (feet, miles, km) in
            % both the world and model system.
            %
            % TO  1204020
            try
                p=all(size(xm)==size(ym));
            catch ME
                error('gridObj/world: size x arg must equal size y arg!');
            end
            if ~p, error('gridObj/world: size x arg must equal size y arg!'); end
                
            xwyw=[xm(:),ym(:)] * [cos(o.anglew) sin(o.anglew); -sin(o.anglew) cos(o.anglew)];

            xw=reshape(xwyw(:,1),size(xm))+o.xw0;
            yw=reshape(xwyw(:,2),size(ym))+o.yw0;
        end
        
        function [xm, ym]=model(o,xw,yw)
            %GRIDOBJ/MODEL -- computes model coordinates form world
            %
            % USAGE:  [xm,ym]=gridObj/model(xw,yw)
            %
            % coordintes xw,yw
            % TO 120420
            try
                p=all(size(xw)==size(yw));
            catch ME
                error('gridObj/model: size xm (first arg) must equal size ym (2nd arg)!');
            end
            if ~p, error('gridObj/model: size xm (first arg) must equal size ym (2nd arg)!'); end
            
            xmym = [xw(:)-o.xw0,yw(:)-o.yw0] * ...
                [cos(o.anglew) -sin(o.anglew); sin(o.anglew) cos(o.anglew)];
            
            xm=reshape(xmym,size(xw));
            ym=reshape(xmym,size(yw));
        end
                
        %% Boundary conditions input using point specifications
        
        function [BCN,PNTSRC]=bcnPoint(o,basename,type,points,vals,conc)
            %GRIDOBJ/BCNPOINT -- generate boundary conditions for points
            %
            % USAGE BCN=gr.bcnPoint(basename,type,points,vals,conc)
            %
            % Points is a list of world coordinates [x y] or [x y z] or {[x
            % y] iLay}. Vals is an [Npoly,n] array with one row per polyling.
            % and as manu values as neceesary, all of who will be appended
            % as extra columns to BCN putputhas. If Vals is a struct array
            % values may also be a string with the name of a header in
            % which will yield a single value for each stress period read from
            % the PER worksheet. Likewise for Conc, where Conc is [Npoly, NCOMP]
            % i.e. with one value for each concentration species involved 
            %
            % TO 120410
            P = o.lineObjects(points);

            zoneArray = gr.const(0);
            zoneArray([P.idx]) = 1;
            
            if ~iscell(vals), vals = {vals}; end
            
            elevation = [P.zm]';
            switch type
                case 'WEL', vals = {elevation vals};
                case 'DRN', vals = {elevation vals};
                case 'GHB', vals = {elevation vals};
                case 'RIV', vals = [elevation vals(1) vals(2) vals(3)];
                otherwise
                    fprintf('%s: unknown TYPE %d, action skipped\n',mfilename,type);
            end

            if nargin<6
                BCN          = bcnZone(basename,type,zoneArray,vals);
            else
                [BCN,PNTSRC] = bcnZone(basename,type,zoneArray,vals,conc);
            end            
        end
        
        %% Boundary conditions input using polyline speicfication
        function [BCN,PNTSRC]=bcnPoly(o,basename,type,poly,vals,conc)
            %GRIDOBJ/BCNPOLY
            %
            % USAGE:
            %    BCN=gr.bcnPoly(basename,type,poly,vals,conc)
            %
            % same as bcnLine all points in the line are inlcuded
            % and vectors of data for a bcnPoly are not allowed, because
            % ambiguous. This may yield huge files, but also effective to
            % put any data into the mode that is available in polyline
            % shapes.
            % TO 130614
            
            zoneArray = double(inpolygon(o.xm,o.ym, poly(:,1),poly(:,2)));
            
            if ~iscell(vals), vals = {vals}; end
            
            switch type
                case 'WEL', vals = [vals(1) vals{2}.*o.AREA(IN)];
                case 'DRN', vals = [vals(1) vals{2}.*o.AREA(IN)];
                case 'GHB', vals = [vals(1) vals{2}.*o.AREA(IN)];
                case 'RIV', vals = [vals(1) vals{2}.*o.AREA(IN) vals(2) vals(3)];
                otherwise
                    fprintf('%s: unknown TYPE %d, action skipped\n',mfilename,type);
            end

            if nargin<6
                BCN          = bcnZone(basename,type,zoneArray,vals);
            else
                [BCN,PNTSRC] = bcnZone(basename,type,zoneArray,vals,conc);
            end            
        end
                
        
        function A=spyBCN(o,BCN)
            %GRIDOBJ/SPYBCN -- spy's BCN where BCN is WEL, CHD, DRN etc
            %
            % USAGE: gr.spyBCN(BN)
            %
            % TO 130101
            
            A=o.const(0);
            if iscell(BCN)
                A(cellIndex(BCN{1}(:,[4,3,2]),o.size))=1;
            else
                A(cellIndex(BCN(:,[4,3,2]),o.size))=1;
            end
            figure; spy(XS(A));
        end
        % Boundary condition input specification through zoneArray
        function [c,h]=streamlines(o,ax,Psi,prange,varargin)
            %GRIDOBJ/STREAMLINES compute streamlines
            %
            % USAGE:  [c,h] = gr.streamlines(ax,Psi,crange,varargin)
            % PSI   = gr.streamlines(ax,Psi); % only computes the PSI with
            %         confining beds
            % plots streamline on ax given Psi, crange and the grid gr
            % TO 120517
            if nargin<3, Psi=ax; end
            
            PSI=NaN(o.Nlay+o.Ncbd+1,size(Psi,2));
            PSI( o.ITlay,:) = Psi(1:end-1,:);
            PSI( o.ITcbd,:) = Psi(find(o.LAYCBD)+1,:);
            PSI(  end   ,:) = Psi( end   ,:);

            if nargin<3
                c = PSI;
            else
                [c,h]=contour(ax,...
                    XS(o.XGR(1,2:end-1,:)),...
                    XS(o.ZGR(1,2:end-1,:)),...
                    PSI,prange,varargin{:});
            end
        end
        function streamlinesUpdate(o,h,Psi)
            %GRIDOBJ/STREAMLINESUPDATE -- update streamlines associated
            %
            % USAGE: gr.streamlinesUpdate(h,Psi)
            %
            % with handle h using Psi as new zdata
            % TO 120517
            PSI=NaN(o.Nlay+o.Ncbd+1,size(Psi,2));
            PSI(o.ITlay,:) = Psi(1:end-1,:);
            PSI(o.ITcbd,:) = Psi(find(o.LAYCBD)+1,:);
            PSI( end   ,:) = Psi( end   ,:);
            set(h,'zdata',PSI);
        end

        function [c,h]=contourXS(o,ax,H,iy,hrange,varargin)
            %GRIDOBJ/CONTOURXS -- contours values in cross section along rows
            %
            % USAGE: gr.contour(ax,H,iy,hrange,varargin)
            %
            % H is assumed to be an array of size(Ny,Nx,Nx) but only the
            % first y plane is contoured.
            %
            % SEE ALSO:  gr.contour gr.contourf gr.contourXS, gr.contourYS gr.contourfXS gr.contourfY
            %
            % TO 120517
            try
                strcmp(get(ax,'type'),'axes');
            catch %#ok<CTCH>
                error('%s: first argument must be an valid axes',mfilename);
            end

           [c,h]=contourf(ax,ones(o.Nz,1)*o.xc,XS(o.ZClay(iy,:,:)),XS(H(iy,:,:)),hrange,varargin{:});
        end

        function [c,h]=contourfXS(o,ax,H,iy,hrange,varargin)
            %GRIDOBJ/CONTOURFXS -- contours values in cross section along rows
            %
            % USAGE: gr.contourf(ax,H,iy,hrange,varargin)
            %
            % H is assumed to be a 3D array of size (Ny,Nx,Nz) but only the
            % first plane is contoured.
            %
            % SEE ALSO:  gr.contour gr.contourf gr.contourXS, gr.contourYS gr.contourfXS gr.contourfY
            %
            % TO 120517
            try
                strcmp(get(ax,'type'),'axes');
            catch %#ok<CTCH>
                error('%s: first argument must be an valid axes',mfilename);
            end
            [c,h]=contourf(ax,XS(o.XC(iy,:,:)),XS(o.ZClay(iy,:,:)),XS(H(iy,:,:)),hrange,varargin{:});
        end

        function [c,h]=contourYS(o,ax,H,ix,hrange,varargin)
            %GRIDOBJ/CONTOURYS -- contours values in cross section along columns
            %
            % USAGE:  gr.contour(ax,H,ix,hrange,varargin)
            %
            % H is assumed to be an array of size(Ny,Nx,Nx) but only the
            % firs x plane is contoured.
            %
            % SEE ALSO:  gr.contour gr.contourf gr.contourXS, gr.contourYS gr.contourfXS gr.contourfY
            %
            % TO 120517
            try
                strcmp(get(ax,'type'),'axes');
            catch %#ok<CTCH>
                error('%s: first argument must be an valid axes',mfilename);
            end

           [c,h]=contour(ax,ones(o.Nz,1)*o.yc,YS(o.ZClay(:,ix,:)),YS(H(:,ix,:)),hrange,varargin{:});
        end

        function [c,h]=contourfYS(o,ax,H,ix,hrange,varargin)
            %GRIDOBJ/CONTOURFYS -- contours values in cross section along columns
            %
            % USAGE: gr.contourf(ax,H,ix,hrange,varargin)
            %
            % H is assumed to be a 3D array of size (Ny,Nx,Nz) but only the
            % first x plane is contoured.
            %
            % SEE ALSO:  gr.contour gr.contourf gr.contourXS, gr.contourYS gr.contourfXS gr.contourfY
            %
            % TO 120517
            try
                strcmp(get(ax,'type'),'axes');
            catch %#ok<CTCH>
                error('%s: first argument must be an valid axes',mfilename);
            end
            [c,h]=contourf(ax,ones(o.Nz,1)*o.yc,YS(o.ZClay(:,ix,:)),YS(H(:,ix,:)),hrange,varargin{:});
        end
        
        function plot(o,alpha,gamma,cl,lw) 
            %GRIDOBJ/PLOT --  should become perspective view on model (not yet implemented)
            %
            % USAGE: gr.plot(alpha,gamma,cl,lw)
            %
            % plots grid in 3D
            %
            % TOFIX:
            %
            % TO 120410

            
            if nargin<5, lw=1; end
            if nargin<4, cl='k'; end
                        
            eObs=[sin(alpha*pi/180),-cos(alpha*pi/180),cos(gamma*pi/180)];
            % nomal vectoren e(x,y,z) to the 6 planes + plane index
            e=[-1  0  0       1    1    1 o.Ny     1 o.Nz;
                1  0  0    o.Nx o.Nx    1 o.Ny     1 o.Nz;
                0 -1  0       1 o.Nx    1    1     1 o.Nz;
                0  1  0       1 o.Nx o.Ny o.Ny     1 o.Nz;
                0  0 -1       1 o.Nx    1 o.Ny  o.Nz o.Nz;
                0  0  1       1 o.Nx    1 o.Ny     1    1];
               
            for i=1:size(e,1)
                if cdot(e(1,1:3),eObs)<0
                   surface(o.XGr(e(i,4:5),e(i,6:7),e(i,8:9)),...
                           o.YGr(e(i,4:5),e(i,6:7),e(i,8:9)),...
                           o.ZGr(e(i,4:5),e(i,6:7),e(i,8:9)),'color',cl,'linwidth',lw);
                end
            end
        end

    end
    methods (Static)
            function warn(Param)
            warning(['gridObj:' Param ':Use' Param 'lay'],...
                'gridObj/%s: LAYCBD~=0, use %slay instead of %s',Param,Param,Param);
            end
    end

end


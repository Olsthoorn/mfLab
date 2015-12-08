classdef gridObj
    %   o=gridObj(xGr,yGr[,zGr [,LAYCBD[, MINDZ[, AXIAL]]])
    %   grid object with methods
    %   Essential is the difference between layers (= model cell layers, LAY)
    %   and confining beds (CBD), not counting as layers but optionally
    %   connect to the bottom of layers. The actual configuration is given
    %   in LAYCBD(Nlay,1) having value >0 when a layer has a CBD below it
    %   and having 0 if not. LAYCBD(end)=0 is guaranteed by mfLab.
    %   Z is the elvation of all tops and bottom of layers and confining beds
    %   combined and, therefore, its third dimension is
    %        Nlay+1<=size(Z,3)<2NLay
    %   To avoid confusion, Nz is hidden. Use Nlay and Ncbd instead.
    %   use zTop and zBot, LAYTOP and LAYBOT to get the elevations of the
    %   tops and bottoms of the model layers and. likewise, use CBDTOP and
    %   CBDBOT to get the tops and bottoms of the confining beds. These,
    %   may however be empty if LAYCBD is all zeros.
    %
    %   TO 110810; gridobject
    
    properties (Constant) % also physically stored
        type='gridObj';
    end
    properties
        MINDZ=0.001;
        AXIAL=0;
        LAYCBD=0;
        isLay; 
        xGr; yGr,                % grid line coordinates, vectors of size Nx+1, Ny+1, Nz+1
        xw0=0; yw0=0; zw0=0; anglew=0;  % position of model zero in real world coordinates
        layersAreUniform,          % true of all layers are uniform
    end
    properties (Access=protected)
        zGrfull
        Zfull
    end
    properties (Access=protected, Dependent=true)
        Nz, dz, zm,
        DX, DY, DZ,
        XM, YM, ZM,
        XGR,YGR,ZGR,
    end
    properties (Dependent=true) % computed when needed bu behavior as if existing
        size, sizeLay, sizeCBD, sizeFULL, % 3rd dimension size resepectively: Nlay Nlay Ncbd Nz
        Nxy, Nxyz              % Nx*Nt,   Nx*Ny*Nlay (total nr of model cells)
        
        ITlay,ITcbd,IBlay,IBcbd, % pointer to top of layer  LAYTOP=Z(:,:,Itop)
        
        Nx, Ny, Nlay, Ncbd,    % size of cells
        
        zLay, zGr,          % derived from zGr
        zTlay,zTcbd,
        zBlay,zBcbd;
        
        xm, ym,zmlay,zmcbd;
        dx, dy
        Xm,Ym,         %2D coordinates of top layer
        
        zMlay, zMcbd,  % center of cells (vectors of size Nx,Ny,Nlay,Ncbd) 
        dzlay, dzcbd,  % size of cells (vectors of size Nx,Ny,Nlay,Ncbd)
        
        XMlay, XMcbd,  % separately for layers and confining beds
        YMlay, YMcbd,  % [Ny,Nx,Nlay] and [Ny,Nx,Ncbd] resp.
        ZMlay, ZMcbd,
        
        DXlay, DXcbd,
        DYlay, DYcbd,
        DZlay, DZcbd,
                
        Z,Zlay,  Zcbd,           % same as ZGRlay and ZGRcbd
        ZGRlay,ZGRcbd,
        ZTlay, ZBlay,
        ZTcbd, ZBcbd,

        Vlay, vlay             % volume of model cells and of entire model
        Vcbd, vcbd             % volume of CBD per cell and total
        AREA, area             % surface area of model cells and entire model
        
        r, rm, dr,             % distance to xm=0 (along the x-axis) for axi-symmetric situations
        xh, yh, zh,            % x to plot heads  [xGr(1) xm(2:end-1) xGr(end)] same for z
        xc, yc, zc,            % same as xc and zc for contouring concentrations
        xp, yp, zp             % [xGr([2:end-1])  zGr all cells + cbds for contouring 
        R, RM                  % distance to center of model
        
        TWOPIR                 % 2*pi*(abs(XM)).^2, full size
end
    methods   
        function o=gridObj(xGr,yGr,Z,LAYCBD,MINDZ,AXIAL)
            if nargin==0; return; end
            if nargin<3,
                error('not enough input arguments for gridObj: use gridObj(xGr,yGr,zGR[,LAYCBD[,MINDX[,AXIAL]]])');
            end
            
            if isempty(xGr), display(xGr); error('gridObj/gridObj call has empty first argument xGr, check it in mf_adapt.'); end
            if isempty(yGr), display(yGr); error('gridObj/gridObj call has empty  2nd  argument yGr, check it in mf_adapt.'); end
            if isempty(Z)  , display(Z  ); error('gridObj/gridObj call has empty  3rd  argument zGr, check it in mf_adapt.'); end
            
            if nargin>3, o.LAYCBD=LAYCBD; end
            if nargin>4, o.MINDZ =MINDZ;  end
            if nargin>5, o.AXIAL =AXIAL;  end
            
            o.xGr = unique(xGr(:))';
            o.yGr = unique(yGr(:));   % flipud is not really necessary
            o.yGr=flipdim(unique(yGr(:)),1);
            zisVector = sum([size(Z,1)==1,size(Z,2)==1,size(Z,3)==1])==2;
            
            if zisVector,
                o.layersAreUniform=1;
                o.Zfull   = XS(flipdim(unique(Z(:)),1));
                o.zGrfull = XS(flipdim(unique(Z(:)),1));
            elseif size(Z,3)==1, % non vector Z must be 3D
                    error('gridObj:gridObj:Znot3D',...
                        '3rd argument Z (or zGr) must be 3D unless it is given as a vector');
            else
                o.zGrfull = mean(mean(Z,1),2);
                o.layersAreUniform = all(Z(1,1,:)==o.zGrfull);
                if o.layersAreUniform
                    o.Zfull=o.zGrfull;
                else
                    o.Zfull=Z;
                end
            end

            %1: guarantee size(..,i) is >=2 in all three i directions
            if numel(o.xGr)==1 || ~isvector(o.xGr),
                error('gridObj/gridObj: arg 1, xGr must be a vector with >1 elements');
            end
            if numel(o.yGr)==1 || ~isvector(o.yGr),
                error('gridObj/gridObj: arg 2, yGr, must be a vector with >1 elmements.');
            end
            if numel(o.zGrfull)==1
                error('gridObj/gridObj: arg 3, (zGr or Z), must have >1 elements');
            end
            
            %3: make sure Z runs from high to low
            % find a point in the grid where Z are not NaN
            Ntrials=10;
            I=find(~isnan(o.Z(:,:,1)),Ntrials,'first');
            for itr=1:Ntrials
                if ~isnan(o.Z(I(itr)+o.Ny*o.Nx*(size(o.Z,3)-1)))
                    break;
                elseif itr==Ntrials
                    error('gridObj:gidObj:zBotNaN',...
                        'In %d trials, I found no z columns with bottom values not NaN',Ntrials);
                end
            end
            RCL = cellIndices(I(itr),[o.Ny o.Nx 1],'RCL');
            if o.Z(RCL(1),RCL(2),end)>o.Z(RCL(1),RCL(2),1),
                o.Z=flipdim(o.Z,3);  % flip Z so that it will run form high to low
            end;

            %4: Guarantee minimum layer thickness
            minLayerThickness=min(min(min(-diff(Z,1,3))));
            if minLayerThickness<o.MINDZ
                warning('gridObj:gridObj:layerThicness',...
                    ['gridObj/gridObj: min layer thickness = %g (<%s),\n',...
                     'correct this before restarting.'],minLayerThickness,o.MINDZ);
            end

            %6: In case size(Z)=[Ny+1,Nx+1,Nz+1] make it [Ny,Nx,Nz]
            if ~o.layersAreUniform
                if size(o.Z,1)>o.Ny,  o.Z=0.5*(o.Z(1:end-1,:,:)+o.Z(2:end,:,:)); end
                if size(o.Z,2)>o.Nx,  o.Z=0.5*(o.Z(:,1:end-1,:)+o.Z(:,2:end,:)); end
            end
            
            if o.AXIAL
                i=find(o.xGr<0,1,'last');
                j=find(o.xGr>0,1,'first');
                if ~isempty(i) && ~isempty(j) && isempty(o.xGr==0) % we must split xGr and Z                    
                    o.xGr=[o.xGr(1:i  )         0            o.xGr(j:end)];
                    o.Z  =[o.Z(:,1:i,:) mean(o.Z(:,i:j,:),2) o.Z(:,j:end,:)];
                end
            end
            
            [o.isLay,o.LAYCBD]=isLayer(o.Nz,o.LAYCBD);
                        
        end
        
        function o = convert(o,Ix,Iy,LAYCBD)
            %% Convert the grid to reflect then new outcut defined by Ix,Iy
            %  where Ix and Iy are the old indices.
            if max(Ix)>o.Nx || min(Ix)<1,
                error('gridObj/convert: Ix indices must be between 1 and %d.',o.Nx); end
            if max(Iy)>o.Ny || min(Iy)<1,
                error('gridObj/convert: Iy indices must be between 1 and %d.',o.Ny); end
            
            if nargin<4, LAYCBD=o.LAYCBD; end
            
            if length(LAYCBD)+sum(LAYCBD>0)~=length(o.LAYCBD)+sum(o.LAYCBD>0),
                error(['gridObj/convert: new and old LAYCBD have different total number of layers\n',...
                     'namely new: Nlay=%d Ncbd=%d, versus old: Nlay=%d Ncbd=%d'],...
                     lenght(LAYCBD),sum(LAYCBD>0),o.Nlay+o.Ncbd);
            end
            
            o.Zfull  = o.Zfull(Iy,Ix,:);
            o.xGr = o.xGr([ Ix(1:end), Ix(end)+1 ]);
            o.yGr = o.yGr([ Iy(1:end); Iy(end)+1 ]);
            [o.isLay o.LAYCBD] = isAquifer(o.Nlay+o.Ncbd,LAYCBD);
        end
        
        function AXIAL = get.AXIAL(o), AXIAL = o.AXIAL;  end
        function LAYCBD= get.LAYCBD(o), LAYCBD=o.LAYCBD; end
        
        function  size    = get.size(o),    size    = [o.Ny, o.Nx, o.Nlay]; end
        function  sizeLay = get.sizeLay(o), sizeLay = [o.Ny, o.Nx, o.Nlay]; end
        function  sizeCBD = get.sizeCBD(o), sizeCBD = [o.Ny, o.Nx, o.Ncbd]; end
        function  sizeFULL= get.sizeFULL(o), sizeFULL=[o.Ny, o.Nx, o.Nz  ]; end
        
        function Nxy  = get.Nxy(o),  Nxy = o.Nx*o.Ny;            end
        function Nxyz = get.Nxyz(o), Nxyz= o.Nx*o.Ny*o.Nlay;     end
        function Nx   = get.Nx(o) ,  Nx   = numel(o.xm);   end
        function Ny   = get.Ny(o) ,  Ny   = numel(o.ym);   end
        function Nz   = get.Nz(o) ,  Nz   = numel(o.zm);   end
        function Nlay = get.Nlay(o), Nlay = sum( o.isLay); end
        function Ncbd = get.Ncbd(o), Ncbd = sum(~o.isLay); end
        
        function xGr   = get.xGr(o), xGr = o.xGr; end
        function yGr   = get.yGr(o), yGr = o.yGr; end
        function zGr   = get.zGr(o), zGr = o.zGrfull; end
        function zLay  = get.zLay(o)
            zLay= cat(3,o.zGr(o.isLay),o.zGr(end));
        end
        function zTlay = get.zTlay(o),zTlay = o.zGr(o.ITlay);   end
        function zBlay = get.zBlay(o),zBlay = o.zGr(o.ITlay+1); end
        function zTcbd = get.zTcbd(o),zTcbd = o.zGr(o.ITcbd);   end
        function zBcbd = get.zBcbd(o),zBcbd = o.zGr(o.ITcbd+1); end
        
        function XGR = get.XGR(o), XGR = repmat(o.xGr,[o.Ny+1,1,o.Nz+1]);  end
        function YGR = get.YGR(o), YGR = repmat(o.yGr,[1,o.Nx+1,o.Nz+1]);  end
        function ZGR = get.ZGR(o), ZGR = o.Zfull; end
        function Zlay= get.Zlay(o),Zlay= cat(3,o.Z(:,:, o.isLay), o.Z(:,:,o.Nz+1)); end
        function Zcbd= get.Zcbd(o),Zcbd= cat(3,o.Z(:,:,~o.isLay), o.Z(:,:,o.Nz+1)); end
        
        function Z    = get.Z(o)
            if o.layersAreUniform
                Z = repmat(o.zGr,[o.Ny,o.Nx]);
            else
                Z   = o.Zfull;
            end 
        end
        function ZTlay = get.ZTlay(o), ZTlay=o.Z(:,:,o.ITlay  ); end 
        function ZBlay = get.ZBlay(o), ZBlay=o.Z(:,:,o.ITlay+1); end 
        function ZTcbd = get.ZTcbd(o), ZTcbd=o.Z(:,:,o.ITcbd  ); end 
        function ZBcbd = get.ZBcbd(o), ZBcbd=o.Z(:,:,o.ITcbd+1); end 
      
        function dx     = get.dx(o),    dx     =  diff(o.xGr,1,2);      end
        function dy     = get.dy(o),    dy     = -diff(o.yGr,1,1);      end
        function dz     = get.dz(o),    dz     = -diff(o.zGr,1,3);      end
        function dzlay  = get.dzlay(o), dzlay  = o.dz(     o.isLay);    end
        function dzcbd  = get.dzcbd(o), dzcbd  = o.DZ(:,:,~o.isLay);    end

        function DX     = get.DX(o),    DX     = repmat(o.dx ,[o.Ny,1,o.Nz]); end
        function DXlay  = get.DXlay(o), DXlay  = o.DX(:,:, o.isLay); end
        function DXcbd  = get.DXcbd(o), DXcbd  = o.DX(:,:,~o.isLay); end
        function DY     = get.DY(o),    DY     = repmat(o.dy ,[1,o.Nx,o.Nz]); end
        function DYlay  = get.DYlay(o), DYlay  = o.DY(:,:, o.isLay); end
        function DYcbd  = get.DYcbd(o), DYcbd  = o.DY(:,:,~o.isLay); end
        function DZ     = get.DZ(o)
            if isvector(o.Z)
                DZ = repmat(-diff(o.Z,1,3),[o.Ny,o.Nx,1]);
            else
                DZ = diff(-o.Z,1,3);
            end
        end
        function DZlay  = get.DZlay(o), DZlay  = o.Z(:,:,o.ITlay)-o.Z(:,:,o.ITlay+1);  end
        function DZcbd  = get.DZcbd(o), DZcbd  = o.Z(:,:,o.ITcbd)-o.Z(:,:,o.ITcbd+1);  end
        
        function xm    = get.xm(o)    , xm    = 0.5*(o.xGr(1:end-1)+o.xGr(2:end));  end
        function ym    = get.ym(o)    , ym    = 0.5*(o.yGr(1:end-1)+o.yGr(2:end));  end
        function zm    = get.zm(o)    , zm    = 0.5*(o.zGr(1:end-1)+o.zGr(2:end));  end  % center of model layers
        function zmlay = get.zmlay(o) , zmlay = o.zm( o.isLay);  end  % center of model layers
        function zmcbd = get.zmcbd(o) , zmcbd = o.zm(~o.isLay);  end  % center of model layers
        
        function zMlay = get.zMlay(o) , zMlay = o.zm( o.isLay);  end  % center of model layers
        function zMcbd = get.zMcbd(o) , zMcbd = o.zm(~o.isLay);  end  % center of model layers

        function Xm    = get.Xm(o),   Xm    = repmat(o.xm,[o.Ny,1]); end
        function Ym    = get.Ym(o),   Ym    = repmat(o.ym,[1,o.Nx]); end
        function XM    = get.XM(o),   XM    = repmat(o.xm ,[o.Ny,1,o.Nz]); end
        function XMlay = get.XMlay(o),XMlay = o.XM(:,:, o.isLay); end  % center Lay+CBD        
        function XMcbd = get.XMcbd(o),XMcbd = o.XM(:,:,~o.isLay); end  % center Lay+CBD        
        function YM    = get.YM(o),   YM    = repmat(o.ym ,[1,o.Nx,o.Nz]); end
        function YMlay = get.YMlay(o),YMlay = o.YM(:,:, o.isLay); end  % center Lay+CBD
        function YMcbd = get.YMcbd(o),YMcbd = o.YM(:,:,~o.isLay); end  % center Lay+CBD        
        function ZM    = get.ZM(o)
            if isvector(o.Z)
                zm_ = 0.5*(o.zGr(1:end-1)+o.zGr(2:end));
                ZM    = repmat(zm_ ,[o.Ny,o.Nx,1]);  % center Lay+CBD
            else
                ZM    = 0.5*(o.Z(:,:,1:end-1)+o.Z(:,:,2:end));
            end
        end
        function ZMlay = get.ZMlay(o),ZMlay = o.ZM(:,:, o.isLay); end  % center Lay+CBD
        function ZMcbd = get.ZMcbd(o),ZMcbd = o.ZM(:,:,~o.isLay); end  % center Lay+CBD
        
        function Vlay = get.Vlay(o),...
                Vlay = o.DX(:,:,o.isLay).* o.DY(:,:,o.isLay).* o.DZ(:,:,o.isLay);
        end
        function vlay = get.vlay(o), vlay = sum(o.Vlay(:));           end
        function Vcbd = get.Vcbd(o); Vcbd = o.DXcbd.* o.DYcbd.* o.DZcbd;  end
        function vcbd = get.vcbd(o); vcbd = sum(o.Vcbd(:));           end
        
        function AREA = get.AREA(o),  AREA = o.dy * o.dx;    end
        function area = get.area(o),  area = sum(o.AREA(:)); end
        
        % axial symmetric, assuming r along x-axis using abs(x)
        function r  = get.r(o),  r  = sqrt(o.xGr.^2); end
        function rm = get.rm(o), rm = abs(o.xm);      end
        function dr = get.dr(o), dr = abs(o.dx);      end
        function R  = get.R( o), R  = sqrt(o.XGR.^2); end
        function RM = get.RM(o), RM = sqrt(o.XM.^2);  end
 
        function TWOPIR = get.TWOPIR(o), TWOPIR=2*pi*o.RM(:,:,o.isLay); end
        
        % world coordinates
        function xw0     = get.xw0(o), xw0=o.xw0; end
        function yw0     = get.yw0(o), yw0=o.yw0; end
        function zw0     = get.zw0(o), zw0=o.zw0; end
        function anglew  = get.anglew(o), anglew=o.anglew; end
        
        % facilitates plotting
        function xh  = get.xh(o) , xh  =  o.xm;  xh([1 end])=   o.xGr([1 end]);    end  % for heads
        function yh  = get.yh(o) , yh  =  o.ym;  yh([1 end])=   o.yGr([1 end]);    end  % for heads in xy plane
        function zh  = get.zh(o) , zh  = XS(o.zm); zh([1 end])=XS(o.zGr([1 end])); end  % for heads in zx plane
        function xc  = get.xc(o) , xc  = o.xm;  xc([1 end])=   o.xGr([1 end]);     end  % for concentrations  
        function yc  = get.yc(o) , yc  = o.ym;  yc([1 end])=   o.yGr([1 end]);     end  % for concentrations  
        function zc  = get.zc(o) , zc  = XS(o.zm); zc([1 end])=XS(o.zGr([1 end])); end  % for concentrations
        function xp  = get.xp(o) , xp= o.xGr(2:end-1); xp([1 end])=o.xGr([1 end]); end  % for Psi
        function yp  = get.yp(o) , yp= o.yGr;                                      end  % for Psi
        function zp  = get.zp(o) , zp= XS(o.zGr);                                  end  % for Psi
        
        function ITlay = get.ITlay(o), ITlay = find( o.isLay); end
        function IBlay = get.IBlay(o), IBlay = o.ITlay+1;      end
        function ITcbd = get.ITcbd(o), ITcbd = find(~o.isLay); end
        function IBcbd = get.IBcbd(o), IBcbd = o.ITcbd+1;      end

        function const = const(o,value)
            % gridObj/const:  array = gr.const(value or vector)
            % if value is scalar:
            % generates a 3D array of size [gr.Ny,gr.Nx,gr.Nlay], with all values=value.
            % if value is vector,
            % generates a 3D array of size [gr.Ny,gr.Nx,length(value)], with
            % array(:,:,i)=value(i);
            %
            % TO 120410
            
            if isscalar(value)
                const=ones(o.Ny,o.Nx,o.Nz)*value;
            else
%                 if numel(value)~=o.Nlay
%                     error('%s/array(..) requires one argument, a scalar or a vector length Nz\n',mfilename);
%                 end
                const=repmat(XS(value(:)),[o.Ny,o.Nx,1]);
            end
        end
        
        function dist = dist(o,x,y)
            % gridObj/dist: dist=gr.dist(x,y)
            % yields the distance of all cell centers to point x,y
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
            % gridObj/setWorld: obj=gr.setWorld(xw0,yw0,anglew)
            % sets world coordinate system for this grid. The
            % coordinates xw0,yw0 match model coordinates 0,0
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
        
        function [xw yw]=world(o,xm,ym)
            % gridObj/world:  [xw,yw]=gr.world(xm,ym)
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
        
        function [xm ym]=model(o,xw,yw)
            % [xm,ym]=gridObj/model(xw,yw) -- computes model coordinates form world
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
        
        function BCN=bcnPoint(o,basename,type,points,vals,Conc)
            % gridObj/bcnPoint: BCN=gr.bcnPoint(basename,type,points,vals,conc)
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
            
            BCN='Stub: should become boundary condition input from points';
        end
        
        %% Boundary conditions input using polyline speicfication
        function BCN=bcnPoly(o,basename,type,poly,vals,conc)
            % gridObj/bcnPoly: BCN=gr.bcnPoly(basename,type,poly,vals,conc)
            % same as bcnLine but now all points in the line are inlcuded
            % and vectors of data for a bcnPoly are not allowed, because
            % ambiguous. This may yield huge files, but also effective to
            % put any data into the mode that is available in polyline
            % shapes.
            BCN='Stub should become boundary condition input from polylines';
        end
        
        %% Boundary condition input using 3D surfaces
        function BCN=bcnSurf(o,basename,type,surf,vals,conc)
            % gridObj/bcnSurf: BCN=gr.bncSurf(basename,surf,vals,conc)
            % not yet implemented
            % TODO:
            BCN='Stub: should become boundary condition input from surfaces';
        end
        
        %% Boundary conditions input using line specificaiton
        function BCN=bcnLine(o,basename,type,line,vals,conc)
            % gridObj/bcnLine: BCN=gr.bcnLine(basename,type,lline,vals,conc)
            % Generate input of boundary conditions given as lines.
            % eg  ...(basename,'WEL',lines,values,conc
            % lines may be a struct with lines x,y,z
            % values must have a value per line of a value for each point
            % of the line or a string with the header of a column in the
            % PER sheet to get values per stress perio. (one value for the
            % whole line for each stress period).
            % lines, as singe array [x y] or [x y z] or as {[x y] iLay}
            % vals as [val1 val1 val3 ...; next line; ...]
            % may also be a {val1 val2 val3 ...; netx line....} 
            % where a val may be a scalar, a vector of length line or
            % a name of a header in the PER sheet in basename
            % coordinates. Facility to input lines and polylines
            % type one of {'WEL','DRN','RIV','GHB','CHD' 'DRT'}
            % not completely developed
            %
            % TO 120413
            
            pline_xyz=line(:,1:3);
            
            P=linegrid(pline_xyz,o.xGr,o.yGr,o.zGr,o.LAYCBD); 
                        
            %% take LAYCBD into account
            for i=1:length(P)
                P(i).iz =layer2aquif(P(i).iz,o.LAYCBD,o.Nz);
                P(i).Idx=(o.Ny*o.Nx)*(P(i).iz-1)+o.Ny*(P(i).ix-1)+P(i).iy;
            end
            
            % s along points
            dS=sqrt(diff(line(:,1)).^2+diff(line(:,2)).^2+diff(line(:,3)).^2);
            S =[0 cumsum(dS)];
            
            % s along the cell centers
            ds=sqrt(diff([P.xm]).^2+diff([P.ym]).^2+diff([P.zm]).^2);
            s=[0 cumsum(ds)];
            
            iPer = 1; % dummy stress period number
            switch lower(type)
                case 'drn'
                    Hd  =interp1(S,line(:,3),s);
                    Lkn =interp1(S,line(:,4),s);
                    Lkn = Lkn.*[P.L];
                    BCN=[ones(size(P))' * iPer; [P.iz]; [P.iy]; [P.ix]; Hd; Lkn ]';
            end

        end
        
        % Boundary condition input specification through zoneArray


        function plotgrid(o)
            % gridObj/plotgrid: gr.plotgrid([clr[,well[,figname[,figcoords]]]])
            % plots the grid
            %
            % TO 120410
            
            plotgrid(o.xGr,o.yGr);
        end
        
        function plot(o,alpha,gamma,cl,lw)  % should become perspective view on model (not yet implemented)
            % gridObj/plot  gr.plot(alpha,gamma,cl,lw)
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
end


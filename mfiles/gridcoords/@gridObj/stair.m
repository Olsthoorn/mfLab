function stair(o,ax,j,Llay,Lcbd,LAYVAR,CBDVAR,varargin)
    % stair(o,ax,j,Llay,Lcbd,LAYVAR,CBDVAR,varargin)
    %   plotting one section along the x-axis, through row j.
    %   used in plotXSec
    %
    %   ax = axis handle
    %   j  = row through which the cross section runs
    %   Llay logical vector indicating the layer to be plotted
    %   Ldbd same for confining beds
    %   LAYVAR a 3D layer variable to be plotted
    %   CBDVAR a 3D confining bed variable to be plotted.
    %   varargin{:} options, see gridObj.
    %
    % SEE ALSO: gridObj.smooth gridObj.plotXSec
    %
    % TO 130321 151207
        
    [hcolor,varargin] = getProp(varargin,'hlines',[]);
    [vcolor,varargin] = getProp(varargin,'vlines',[]);
    
    x           = NaN(1,2*o.Nx);
    x(1:2:end-1)= o.xGr(1:end-1);
    x(2:2:end)  = o.xGr(2:end);

    zbL = NaN(o.Nlay,2*o.Nx);
    zbL(:,1:2:end-1) = XS(o.ZBlay(j,:,:));
    zbL(:,2:2:end  ) = XS(o.ZBlay(j,:,:));

    ztL = NaN(o.Nlay,2*o.Nx);
    ztL(:,1:2:end-1) = XS(o.ZTlay(j,:,:));
    ztL(:,2:2:end  ) = XS(o.ZTlay(j,:,:)); 
                           
    x2 = [x x(end:-1:1)];
    
    %% first plot patches
    if any(Llay)
        if ischar(LAYVAR)
            for iLay=find(Llay(:)')
                z2 = [zbL(iLay,:) ztL(iLay,end:-1:1)];
                color = LAYVAR(min(iLay,numel(LAYVAR)));
                fill(x2,z2,color(1),'parent',ax,varargin{:});
            end
        else
            C = NaN(o.Nlay,numel(x));
            C(:,1:2:end-1) = XS(LAYVAR(j,:,:));
            C(:,2:2:end  ) = XS(LAYVAR(j,:,:));
            for iLay=find(Llay(:)')
                z2 = [zbL(iLay,:) ztL(iLay,end:-1:1)];
                fill(x2,z2,[C(iLay,:) C(iLay,end:-1:1)],'parent',ax,varargin{:});
            end
        end
    end    
    
    if any(Lcbd)
        
        ztC = NaN(o.Ncbd,numel(x));
        zbC = NaN(o.Ncbd,numel(x));
        
        ztC(:,1:2:end-1) = XS(o.ZTcbd(j,:,:));
        ztC(:,2:2:end  ) = XS(o.ZTcbd(j,:,:));
        zbC(:,1:2:end-1) = XS(o.ZBcbd(j,:,:));
        zbC(:,2:2:end  ) = XS(o.ZBcbd(j,:,:));
        
        if ischar(CBDVAR)
            for iCbd=find(Lcbd(:)')
                color = CBDVAR(min(iCbd,numel(CBDVAR)));
                fill(x2,[zbC(iCbd,:) ztC(iCbd,end:-1:1)],color,'parent',ax,varargin{:});
            end
        else
            C = NaN(o.Ncbd,numel(x));
            C(:,1:2:end-1) = XS(CBDVAR(j,:,:));
            C(:,2:2:end  ) = XS(CBDVAR(j,:,:));
            for iCbd=find(Lcbd(:)')
                fill(x2,[zbC(iCbd,:) ztC(iCbd,end:-1:1)],[C(iCbd,:) C(iCbd,end:-1:1)],'parent',ax,varargin{:});
            end
        end
    end
    
   %% plot lines last over patches
   if ~isempty(hcolor)
        for i=size(zbL,1):-1:1
            plot(ax,x,zbL(i,:),'color',hcolor);
        end
    end
    
    if ~isempty(vcolor)
        for ix=1:numel(x)
            plot(ax,x([ix ix]),[zbL(end,ix) ztL(1,ix)],'color',vcolor);
        end
    end
end


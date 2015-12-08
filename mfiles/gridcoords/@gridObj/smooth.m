function smooth(o,ax,jRow,Llay,Lcbd,LAYVAR,CBDVAR,varargin)
    % smooth(o,ax,jRow,Llay,Lcbd,LAYVAR,CBDVAR,varargin)
    %   plotting one section along the x-axis, through row jRow.
    %   used in plotXSec
    %
    %   ax   = axis handle
    %   jRow = row through which the cross section runs
    %   Llay logical vector indicating the layer to be plotted
    %   Icbd same for confining beds
    %   LAYVAR a 3D layer variable to be plotted
    %   CBDVAR a 3D confining bed variable to be plotted.
    %   varargin{:} options, see gridObj.
    %
    % SEE ALSO: gridObj.stair gridObj.plotXSec
    %
    % TO 130321 151207
        
    grey = [0.8 0.8 0.8];
    
    zbL = XS(o.ZBlay(jRow,:,:));
    ztL = XS(o.ZTlay(jRow,:,:));

    % look for hlines in varargin
    [hcolor, varargin] = getProp(varargin,'hlines',[]);
    [vcolor, varargin] = getProp(varargin,'vlines',[]);
    
    if ischar(LAYVAR)
       for iLay=find(Llay(:)')
          color = LAYVAR(min(iLay,numel(LAYVAR))); 
          fill([o.xc o.xc(end:-1:1)],[zbL(iLay,:) ztL(iLay,end:-1:1)],color,'parent',ax,varargin{:});
       end
    else
        for iLay=find(Llay(:)')
            color = XS(LAYVAR(jRow,:,iLay));            
            fill([o.xc o.xc(end:-1:1)],[zbL(iLay,:) ztL(iLay,end:-1:1)],[color color(end:-1:1)],'parent',ax,varargin{:});
        end
    end
    
    if ~isempty(CBDVAR)
        zbC = XS(o.ZBcbd(jRow,:,:));
        ztC = XS(o.ZTcbd(jRow,:,:));

        if ischar(CBDVAR)        
           for iCbd=find(Lcbd(:)')
              color = CBDVAR(min(iCbd,numel(CBDVAR))); 
              fill([o.xc o.xc(end:-1:1)],[zbC(iCbd,:) ztC(iCbd,end:-1:1)],color,'parent',ax,varargin{:});
           end
        else
            for iCbd=find(Lcbd(:)')
                color = XS(CBDVAR(jRow,:,iCbd));
                fill([o.xc o.xc(end:-1:1)],[zbC(iCbd,:) ztC(iCbd,end:-1:1)],[color color(end:-1:1)],'parent',ax,varargin{:});
            end
        end
    end
    
    %% Lines are plotted last on top of possible patches
    % always plot top and bottom
    plot(ax,o.xc,ztL(  1,:),'color',grey);
    plot(ax,o.xc,zbL(end,:),'color',grey);


    if ~isempty(hcolor)
        plot(ax,o.xc,ztL(1,:),'color',hcolor);
        for iLay = 1:o.Nlay
            plot(ax,o.xc,zbL(iLay,:),'color',hcolor);
        end            
    end
    
    if ~isempty(vcolor)
        x = [o.xGr; o.xGr];
        z = [ztL(  1,1) 0.5*(ztL(  1,1:end-1)+ztL(  1,2:end)) ztL(  1,end); ...
             zbL(end,1) 0.5*(zbL(end,1:end-1)+zbL(end,2:end)) zbL(end,end)];
         
        for ix=1:o.Nx+1
            plot(ax,x(:,ix),z(:,ix),'color',vcolor);
        end
    end
        
 

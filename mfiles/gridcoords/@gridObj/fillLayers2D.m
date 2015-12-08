function h=fillLayers2D(o,K,iRow,ILay,varargin)
    % h = fillLayers2D(ax,gr,K,iRow,ILay,varargin)
    %
    % fill the layers (probably aquitards) in a cross section with colors that reflect their
    % resistance
    %
    % EXAMPLE
    %   h = fillLayers2D(gr,K,[2 5 6],'color','white','alpha',0.3,'parent',ax)
    %
    % TO 120507
    h = NaN(length(ILay),1);

    for iL=1:length(ILay)
      h(iL)= fill([o.xc,                 o.xc(end:-1:1)               ],...
              [o.ZBlay(iRow,:,iL), o.ZTlay(iRow,end:-1:1,iL)],...
              [K(iRow,:,iL) K(iRow,end:-1:1,iL)],varargin{:});
    end
end

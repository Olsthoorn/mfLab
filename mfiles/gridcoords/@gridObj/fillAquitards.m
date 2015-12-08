function h = fillAquitards(o,VKCB,iRow,varargin)
    % h = fillAquitards(ax,gr,VKCB,iRow,varargin)
    %
    % fill the aquitards in a cross section with colors that reflect their
    % resistance
    %
    % EXAMPLE
    %   h = fillAquitards(gr,VKCB,'color','white','alpha',0.3,'parent',ax)
    %
    % TO 120507

    h=NaN(o.Ncbd,1);

    for icbd=1:o.Ncbd
      h(icbd)= fill([o.xc,                 o.xc(end:-1:1)               ],...
              [o.ZBcbd(iRow,:,icbd), o.ZTcbd(iRow,end:-1:1,icbd)],...
              [VKCB(iRow,:,icbd) VKCB(iRow,end:-1:1,icbd)],varargin{:});
    end
end

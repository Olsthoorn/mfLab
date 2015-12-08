function IN = menyanthesIn(varargin)
    %MENYANTHESIN -- plots and gives the heads at pointObj
    %locations, considering them as piezometers
    %
    % USAGE: IN = menyanthesIn(['name',names,]['t0',t0,]['tp'|'prec' prec,]['te'|'evap' evap,]['tep' tep,]['x',xcoord,]['y',ycoord])
    %
    %  names is {'name1','name2','name3','name4',...} etc.
    %  prec  is [t P]
    %  evap  is [t E]
    %  evap  is [t E P]
    %
    % TO 131011
    
    template.ID        =  [];
    template.name       = '';
    template.type       = '';
    template.values     = [];
    template.xcoord     = [];
    template.ycoord     = [];
    template.surflev   =  [];
    template.filtnr     = [];
    template.upfiltlev  = [];
    template.lowfiltlev = [];
    template.datlog_serial  = [];
    
    [t0,    varargin] = getType(varargin,'t0',0);
    [names, varargin] = getProp(varargin,'name',{'name'});    
    [xcoord,varargin] = getProp(varargin,{'x','xcoord'},0);
    [ycoord,varargin] = getProp(varargin,{'y','ycoord'},0);

    [tpe,varargin] = getProp(varargin,'tpe',[]);
    [te, varargin] = getProp(varargin,{'te','evap'},[]);
    [tp,    ~    ] = getProp(varargin,{'tp','prec'},[]);

    k=0;
    if ~isempty(te)
        te(:,1) = te(:,1) + t0;
        k=k+1;
        IN(k)         = template;
        IN(k).ID      = k;
        IN(k).name    = sprintf('EVAP%02d',k);
        IN(k).type    = 'EVAP';
        IN(k).values  = te;
        IN(k).xcoord  = xcoord(min(1,numel(xcoord)));
        IN(k).ycoord  = ycoord(min(1,numel(ycoord)));
    end
    
    if ~isempty(tp)
        tp(:,1) = tp(:,1) + t0;
        k=k+1;
        IN(k)         = template;
        IN(k).ID      = k;
        IN(k).name    = sprintf('PREC%02d',k);
        IN(k).type    = 'PREC';
        IN(k).values  = tp;
        IN(k).xcoord  = xcoord(min(1,numel(xcoord)));
        IN(k).ycoord  = ycoord(min(1,numel(ycoord)));
    end
    
    if ~isempty(tpe)
        tpe(:,1) = tpe(:,1) + t0;        
        k=k+1;
        IN(k)        = template;
        IN(k).ID     = k;
        IN(k).name   = sprintf('PREC%02d',k);
        IN(k).type   = 'PREC';
        IN(k).values = tpe(:,[1 2]);
        IN(k).xcoord  = xcoord(min(1,numel(xcoord)));
        IN(k).ycoord  = ycoord(min(1,numel(ycoord)));
        
        k=k+1;
        IN(k)        = template;
        IN(k).ID     = k;
        IN(k).name   = sprintf('EVAP%02d',k);
        IN(k).type   = 'EVAP';
        IN(k).values = tpe(:,[1 3]);
        IN(k).xcoord  = xcoord(min(1,numel(xcoord)));
        IN(k).ycoord  = ycoord(min(1,numel(ycoord)));

    end
    
end

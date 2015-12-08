function plotXSec(o,varargin)
% h=gr.plotXSec(varargin) -- plots XSections along xaxis,
%   through all Jrow of the mesh. These cross secitons are plotted on their
%   own figure/axis or subplot axis.
%
% Aquitards are filled with a sequence of colors. This is automatically the
% case of confining beds are present, that is, any gr,LAYCBD ~=0.
%
% Aquifers may also be filled if ~isempty(ILay).
%
% 
% USAGE
%  [hlay hcbd] = gr.plotXSec([ax],Jrow,Ilay,options)
%  [hlay hcbd] = gr.plotXSec([ax],Jrow,Ilay,Icbd,options)
%  [hlay hcbd] = gr.plotXSec([ax],Jrow,'all',plotoptions)
%  [hlay hcbd] = gr.plotXSec([ax],Jrow,'lay',Ilay,options)
%  [hlay hcbd] = gr.plotXSec([ax],Jrow,'lay',LAYVAR,'cbd',CBDVAR,'clim',clim,options)
%
% Jrow: (first numeric input)  vector of the cross sections and or layers to be plotted.
% Ilay: (second numeric input) indices or logical vector telling which model layers to plot
% Icbd: (third numeric input)  indices or logical vector telling which confining beds to plot 
%
% options are propertyname property pairs:
%   'figure',figname  --> new figure requested
%   'posiion',position --> postion of figure on screen
%   'hlines',colorSpec  --> plot lines between layer in specified color
%   'vlines',colorSpec  --> plots vertical lines between cells in specified color
%   'lay',3Darray  --> to use for coloring XSec
%   'cbd',3Darray  --> to use for coloring XSec
%
% one-word options to be used are the following:
%   'all'     --> color all layers and confining beds
%   'smooth' --> plot layer interfaces as contniuous lines
%   'stair'  --> plot layer interfaces as stair lines (showing actual cell bottoms and tops)
% plotOptions
%
% plotOptions (any valid options for Matlab's fill and patch functions)
%   'edgecolor',colorSpec
%   'faceColor',colorSpec
%   'edgeAlpha',fraction
%   'faceAlpha',fraction
%   'lineWidth',value
%   'color',colorSpec
%   plus any options that are accepted by matlab's fill and patch
%
% NOTE
% Ilay and Icbd may be vectors of layer indices or logical vectors, i.e. the
% result of a logical expression. The latter allows parameterizing the layers
% to be filled. For instance, by demanding that its conductivity is less than
% some prescribed value.
% For instance a criterion like this
%    max(max(HK(:,:,:),[],1),[],2)<0.03
% is such a logical vector selecting only the layers whose max HK<0.03.
%
%
% EXAMPLE:
%  h = gr.plotXSec(1,[5 7 9],'lay',HK,'cbd',VKCB,'hlines','y','vlines','b','linewidth',2,'facecolor','g');
%  h = gr.plotXSec(1,[5 7 9],'smooth','hlines','y','vlines','b','linewidth',2,'facecolor','g');
%  h = gr.plotXSec(1,[5 7 9],'lay',log10(HK),'hlines',k,'faceAlpha',0.25);
%  h = gr.plotXSec(1,'all','figure','vlines''g','smooth');
%  h = gr.plotXSec(1,mean(mean(HK,2),1)<15,'stair')
%
%     mean(mean(HK,1),2)<15 says that the average layer conductivity<15
%
% SEE ALSO: gr.plotYSec | gr.plotGrid | gr.mesh | gr.plotXS | gr.plotYS
%
% TO 120501 120531 130208 151207
%
% Copyright 2009-2013 Theo Olsthoorn
% under free software foundation GNU license version 3 or later

%%

vararg = {};

fsz = 12;

%% fig name and position
[ax     ,varargin] = getProp(varargin,'axis',[]);
[ax     ,varargin] = getType(varargin,'axis',ax);
[figPos ,varargin] = getProp(varargin,'figpos',screenPos(0.75));
[figPos ,varargin] = getProp(varargin,'pos',figPos);
[figName,varargin] = getProp(varargin,'fig','noName');
[fsz    ,varargin] = getProp(varargin,'fonts',fsz);

newFig = ~strcmp(figName,'noName') && isempty(ax);


%% which type of plot?

[  ~     ,varargin] = getWord(varargin,'smooth');
[isStair ,varargin] = getWord(varargin,'stair');

if isStair
    plotfun = @stair;
else
    plotfun = @smooth;
end

[dir,     varargin] = getProp(varargin,'dir','x');
if ~ischar(dir) || ~ismember(lower(dir),{'x' 'y'})
    error('%s: dir must be ''x'' or ''y''',mfilename);
end


%% plot hlines and or vlines?

[lines ,  varargin] = getProp(varargin,'lines','off');
[hlines,  varargin] = getProp(varargin,'hlines',lines);
[vlines,  varargin] = getProp(varargin,'vlines',lines);

if strcmpi(hlines,'on'), hlines = 'k'; end
if strcmpi(vlines,'on'), vlines = 'k'; end

if ~strcmpi(hlines,'off')
    vararg = [vararg, {'hlines',hlines}];
end
if ~strcmpi(vlines,'off')
    vararg = [vararg, {'vlines',vlines}];
end


%% See which layers have to be plotted

[clim,varargin] = getProp(varargin,'clim',[]);
if ~isempty(clim)
    vararg = [varargin, {'clim',clim}];
end


%% See if specific variables (3D arrays have to be plotted)

LayError = false;

[LAYVAR,varargin] = getProp(varargin,'lay',[]);
if ~isnumeric(LAYVAR), LayError=true; end

[CBDVAR,varargin] = getProp(varargin,'cbd',[]);
if ~isnumeric(CBDVAR), LayError = true; end

if LayError
    error('%s: illegal input, argument after ''lay'' or ''cbd'' must be a numeric 3D array or 3D vector',mfilename);
end


%% title

[ttl,varargin ] = getProp(varargin,'title',mfilename);

%% xlabel

if dir=='x'
    if o.AXIAL        
        [XLABEL,varargin] = getProp(varargin,'xlabel',sprintf(' %s [m] ','r'));
    else
        [XLABEL,varargin] = getProp(varargin,'xlabel',sprintf(' %s [m] ',dir));
    end
    [ZLABEL,varargin] = getProp(varargin,'ylabel',[]);
    if isempty(ZLABEL)
        [ZLABEL,varargin] = getProp(varargin,'zlabel','z [m]');
    end
else
    [XLABEL,varargin] = getProp(varargin,'xlabel','');
    if isempty(XLABEL)
        [XLABEL,varargin] = getProp(varargin,'ylabel',sprintf('SOUTH %s [m] NORTH',dir));
    end
    [ZLABEL,varargin] = getProp(varargin,{'ylabel','zlabel'},'z [m]');
end

%% choose Jrow = 1 if unspecified

[Jrow,varargin] = getNext(varargin,'double',1);

% see if we have to use a new figure?
if numel(Jrow)>1
    newFig =true;
end


%% Llay Lcbd
[allLayers ,varargin] = getWord(varargin,'all');
if allLayers
    Llay = o.ITlay;
    Lcbd = o.ITcbd;
else
    [Llay,varargin] = getNext(varargin,{'logical','double'},o.ITlay);
    [Lcbd,varargin] = getNext(varargin,{'logical','double'},o.ITcbd);
end

%% plot as many cross sections as desired

for jRow=Jrow(Jrow>0 & Jrow<=o.Ny) % for all desired cross sections

    if newFig
        figure('name',[figName sprintf(' %d',jRow)],'position',figPos);
        
        ax = axes('nextplot','add','xgrid','on','ygrid','on','fontsize',fsz);
        if ~isempty(clim)
            set(ax,'clim',clim);
        end
        xlabel(ax,XLABEL,'fontsize',fsz); ylabel(ax,ZLABEL,'fontsize',fsz);
    else
        if isempty(ax), ax=gca; end
        set(ax,'nextplot','add');
    end

    if dir=='x'
        title(ax,sprintf('%s along x-axis at ym = %6g m, row#=%d',ttl,o.ym(jRow),jRow),'fontsize',fsz);
        %fprintf('%s ,section along x-axis: ym = %6g, row=%d\n'       ,ttl,o.ym(jRow),jRow);
        xlabel(ax,XLABEL,'fontsize',fsz);
    elseif  dir=='y'
        title(ax,sprintf('%s along y-axis at xm = %6g m, col#=%d',ttl,o.ym(jRow),jRow),'fontsize',fsz);
        %fprintf('%s ,section along y-axis: xm = %6g, column=%d\n'    ,ttl,o.ym(jRow),o.Ny-jRow+1,'fontsize',fsz);
        xlabel(ax,XLABEL,'fontsize',fsz);
    else
        error('%s: dir must be ''x'' or ''y'' ');
    end

    %% Plotting or coloring of layers if requested
    if isempty(LAYVAR)
        LAYVAR = mf_color(1:o.Nlay);
    elseif ischar(LAYVAR)
        % skip
    elseif isnumeric(LAYVAR)
        if size(LAYVAR,1) ~= o.Ny || size(LAYVAR,2) ~= o.Nx || size(LAYVAR,3) ~= o.Nlay
            error('%s: CBDVAR must be of size [Ny Nx Ncbd] = [%d %d %d], not [%d %d %d].',...
                mfilename,o.size(),size(LAYVAR));
        end
        % skip
    else
        error('%s: LAYVAR must be of type char or numeric not <<%s>>',mfilename,class(LAYVAR));
    end

    %% Plotting confining beds  if requested
    if isempty(CBDVAR)
        CBDVAR = mf_color(1:o.Ncbd)';
    elseif ischar(CBDVAR)
        % skip
    elseif isnumeric(CBDVAR)
        if size(CBDVAR,1) ~= o.Ny || size(CBDVAR,2) ~= o.Nx || size(CBDVAR,3) ~= o.Ncbd
            error('%s: CBDVAR must be of size [Ny Nx Ncbd] = [%d %d %d], not [%d %d %d].',...
                mfilename,o.Ny,o.Nx,o.Ncbd,size(CBDVAR));
        end
    else
        error('%s: CBDVAR must be of type char or numeric not <<%s>>',mfilename,class(CBDVAR));
    end

    plotfun(o,ax,jRow,Llay,Lcbd,LAYVAR,CBDVAR,vararg{:},varargin{:});
end

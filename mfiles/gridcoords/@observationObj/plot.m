function o=plot(o,varargin)
%OBSERVATIONOBJ/PLOT --- plot observation series held by observatinObj
%
% USAGE:
%     first generate observatin objects, see
%     help observationObj
%     the plot  what is desired
%     obs = obs.plot('nameOfParameter',varargin)
%     see obs.UserData for names of stored parameters
%
% Concentrations will be stored as
%     obs.UserData.(speciesName)
%     the weighed value can be retrieved by the list
%     [obs(io).UserData.(speciesName).value]
%
% sorbed concentrations will be in obs.CS(iComp,1:NT);
%
% SEE ALSO: wellObj.setCout
%
% TO 130403

[ax     ,varargin] = getNext(varargin,'axis',[]);
[ax     ,varargin] = getProp(varargin,'axis',ax);
[figName,varargin] = getProp(varargin,'fig','');
[figPos ,varargin] = getProp(varargin,{'pos','figPos'},screenPos(0.6));
[xLbl   ,varargin] = getProp(varargin,'xlabel','time [d]');
[yLbl   ,varargin] = getProp(varargin,'ylabel',' ? ');

% field is the name of the variable to be plotted (H, DDN, name of species)
% see data in observatinObj.UserData
[field  ,varargin] = getProp(varargin,'field',[]);
if isempty(field)
    [field  ,varargin] = getNext(varargin,'char',[]);
end

if isempty(field)
    error('%s: A field to be plotted must be specified. See observation.UserData for fields.',...
        mfilename);
else
    if strfind('heads',lower(field))
        field = 'H';
    elseif strfind('drawdown',lower(field))
        field = 'DDN';
    end
end

%% New figure explicitly requested
if ~isempty(figName)
    figure('name',figName,'position',figPos);
    ax = axes('nextplot','add');
else
    if isempty(get(0,'children'))
        figure('name',figName,'position',figPos);
        ax = axes('nextplot','add');
    else
        set(gcf,'position',figPos);
        if isempty(ax)
            ax = axes('nextplot','add');
        end
    end
end

set(ax,'xgrid','on','ygrid','on');

xlabel(ax,xLbl);
ylabel(ax,yLbl);
title(field);
grid on;

leg{numel(o),1} = 'allocate';
h(numel(o))=NaN;
for io=1:numel(o)
    if isfield(o(io).UserData,field)
        data = o(io).UserData.(field);
        t    = o(io).t;
        h(io) =plot(t,[data.value],o(io).lineSpec,'lineWidth',o(io).lineWidth,varargin{:});
        leg{io} = o(io).legend;
    else
        error('%s: observationObj.UserData has no field <<%s>>',mfilename,field);
    end
end

legend(leg,'location','northeastOutside');


%% eventually pass plot options to gca
if ~isempty(varargin) && ischar(varargin{1})
    for i=1:2:length(varargin)
        try
            set(ax,varargin{i:i+1});
        catch ME
            fprintf('%s\n',ME.message);
            % skip
        end
    end
end



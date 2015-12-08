function well = cOut(o,well,STCONC,varargin)
    % well = animate.cOut(well [,STCONC],varargin); % shows and add cOut to wells
    % if STCONC is given it add Cout-STCONC.
    % Plots this output concentration of wells versus time
    % TO 121016
    
    % wells are required for cOut
    if nargin<2 || (~strcmpi(class(well),'wellObj') && ~strcmpi(class(well),'MNW1Obj'))
        error('%s: this method requires well of class wellObj as second argument.\n',mfilename);
    end
    
    % Use conc or the diff with STCONC ?
    Delta = nargin>2 && ~isempty(STCONC);
        
    if Delta, prefix = '\Delta'; else prefix = ''; end
    
    % varargin may be used to set plotting properties
    if nargin<4 || isempty(varargin), varargin={'visible','on'}; end
    
    %% subtract STCONC if Delta==true
    for iComp=o.NCOMP:-1:1
        C{iComp} = o.UCN{iComp};
        if Delta
            for it=1:numel(C{iComp})
                C{iComp}(it).values = o.UCN{iComp}(it).values - STCONC{iComp};
            end
        end
        % set well.Cout
        well = well.setCout(C{iComp},iComp);
    end
    
    %% plot the concentrations or their difference with STCONC
    for iComp=1:o.NCOMP
                    
        figure('Name',['Output ', prefix, o.titleStr{iComp}, ' of the wells']);
        
        axes('nextplot','add','xgrid','on','ygrid','on');
        
        %% use days or years as time axes ?
        if well(end).t(end)>3650,
            xlabel('time [y]');
            tfac = 365.24;
        else
            xlabel('time [d]');
            tfac = 1;
        end
        
        %% title
        ylabel(o.titleStr{iComp});
        title(['Output ', prefix, o.titleStr{iComp} ' of the wells']);

        leg{numel(well)} = '_'; %#ok
        
        %% plot of all wells, one subplot per species, using color, lineWidth and lineStyle
        % obtained form well(iw).UserData (see table wells in spreadsheet) and or varargin
        for iw=1:numel(well)            
            leg{iw}=sprintf('wel%02d %s  scr:%5.0f %5.0f',...
                well(iw).nr,...
                well(iw).name,...
                well(iw).z([1 end])); %#ok
                plot(well(iw).t(1+numel(well(iw).t)-numel(well(iw).Dt):end)/tfac,well(iw).Cout(iComp,:),...
                            pick(iw,'rbkgmc'),...
                'lineWidth',pick(iw,[1 1.5 2 2.5]),...
                'lineStyle',pick(iw,{'-','--',':','-.'}),...
                varargin{:});
        end
        legend(leg,4);
    end
end

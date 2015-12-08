        function slices(o,gr,it)
            
            if nargin<2, error('%s: needs grid as argument'); end
        
            time = [o.C.time];

            %% make movie
            if nargin<3,
                it = length(time);
            end

            %% Temperature
            figure('Name','Temperature');
            axes('nextplot','add','clim',o.Tlim);
            xlabel('x [m]'); ylabel('elevation [m]');
            title([o.titleStr ', temp, time = ' datestr(time(it),'mmm-yy')]);
            
            for iw=1:length(o.well)
                if o.well(iw).UserData.isCold,
                    o.well(iw).UserData.color = 'b';
                else
                    o.well(iw).UserData.color = 'r';
                end
            end

            for iw=length(o.well):-1:1
                iy = hit(gr.yGr,o.well(iw).y);            
                surf(squeeze(gr.XGR(iy,:,:)),squeeze(gr.YGR(iy,:,:)),squeeze(gr.ZGR(iy,:,:)),squeeze(o.T(it).values(iy,:,:)),'edgecolor','none');
                o.well(iw)=gr.plotWells(gca,o.well(iw));
            end
            colorbar;
            view(3);

            %% Salinity
            figure('Name','Salinity');
            axes('nextplot','add','clim',o.cLim);
            xlabel('x [m]'); ylabel('elevation [m]');
            title([o.titleStr,', conc, t = ' datestr(time(it),'mmm-yy')]);

            for iw=length(o.well):-1:1
                iy = hit(gr.yGr,o.well(iw).y);            
                surf(squeeze(gr.XGR(iy,:,:)),squeeze(gr.YGR(iy,:,:)),squeeze(gr.ZGR(iy,:,:)),squeeze(o.C(it).values(iy,:,:)),'edgecolor','none');
                o.well(iw)=gr.plotWells(gca,o.well(iw));
            end
            colorbar;
            view(3);


            %% Head
            figure('Name','Point-water head');
            axes('nextplot','add','clim',o.hLim);
            xlabel('x [m]'); ylabel('elevation [m]');
            title([o.titleStr ', PW heads, t = ',datestr(time(it),'mmm-yy')]);

            for iw=length(o.well):-1:1
                iy = hit(gr.yGr,o.well(iw).y);            
                surf(squeeze(gr.XGR(iy,:,:)),squeeze(gr.YGR(iy,:,:)),squeeze(gr.ZGR(iy,:,:)),squeeze(o.H(it).values(iy,:,:)),'edgecolor','none');
                o.well(iw)=gr.plotWells(gca,o.well(iw));
            end
            colorbar;
            view(3);
        end

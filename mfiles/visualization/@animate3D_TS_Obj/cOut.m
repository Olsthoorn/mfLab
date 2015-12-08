        function cOut(o,debugLevel)
            % Plot output concentration of wells versus time
            o.well = o.well.setCout(o.C,2);
            o.well = o.well.setCout(o.T,1);

            if nargin>1,
                dbtxt = sprintf(', debugLevel = %d',debugLevel);
            else
                dbtxt ='';
            end
            figure('Name',['Output temp and concentration of wells',dbtxt]);
            
            o.ax(1)=subplot(2,1,1,'nextplot','add','xgrid','on','ygrid','on');
            o.ax(2)=subplot(2,1,2,'nextplot','add','xgrid','on','ygrid','on');

            ylabel(o.ax(1),'temperature');
            ylabel(o.ax(2),'concentration');

            title(o.ax(1),[o.titleStr 'temperature output of wells']);
            title(o.ax(2),[o.titleStr 'concentration output of wells']);
            
            xlabel(o.ax(1),[datestr(o.well(1).t(1),'yyyy') '-' datestr(o.well(1).t(end),'yyyy')]);
            xlabel(o.ax(2),[datestr(o.well(1).t(1),'yyyy') '-' datestr(o.well(1).t(end),'yyyy')]);
            
            leg{length(o.well)} = '';
            for iw=1:length(o.well)
                leg{iw}=sprintf('well%02d scr %.0f %.0f',o.well(iw).nr,o.well(iw).z);
                linewidth = o.well(iw).UserData.linewidth;
                linecolor = o.well(iw).UserData.color;
%                if o.well(iw).UserData.type =='X', linestyle='--'; else linestyle ='-'; end % dummy
                plot(o.ax(1),o.well(iw).t(2:end),o.well(iw).Cout(1,:),linecolor,'linewidth',linewidth);
                plot(o.ax(2),o.well(iw).t(2:end),o.well(iw).Cout(2,:),linecolor,'linewidth',linewidth);
            end
            datetick(o.ax(1));
            datetick(o.ax(2));
            
            legend(o.ax(1),leg);
            legend(o.ax(2),leg);
        end

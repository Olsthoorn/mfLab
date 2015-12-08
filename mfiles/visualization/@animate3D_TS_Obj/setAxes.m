        function o=setAxes(o,gr,setlim)
            
            if exist('setlim','var') && ~isempty(setlim)
                L=200;
                o.xLim = [min([o.well.x])-L, max([o.well.x])+L];
                o.yLim = [min([o.well.y])-L, max([o.well.y])+L];
            else
                o.xLim = gr.xGr([1 end]);
                o.yLim = gr.yGr([end 1]);

            end

            %% Axes positions            
            for ia = 6:-1:1
                o.ax(ia) = subplot(3,2,ia ,'xgrid','on','ygrid','on','zgrid','on','nextplot','add');
                colorbar
                xlabel(o.ax(ia),'x [m]');
                ylabel(o.ax(ia),'y [m]');
                zlabel(o.ax(ia),'z [m]');
                
                 switch ia
                     case 1, set(o.ax(ia),'clim',o.trange([1 end]),'xlim',o.xLim,'ylim',o.yLim);
                     case 2, set(o.ax(ia),'clim',o.crange([1 end]),'xlim',o.xLim,'ylim',o.yLim);
                     case 3, set(o.ax(ia),'clim',o.trange([1 end]),'xlim',o.xLim);
                     case 4, set(o.ax(ia),'clim',o.crange([1 end]),'xlim',o.xLim);
                     case 5, set(o.ax(ia),'clim',o.trange([1 end]),'xlim',o.xLim,'ylim',o.yLim); 
                     case 6, set(o.ax(ia),'clim',o.crange([1 end]),'xlim',o.xLim,'ylim',o.yLim);
                 end
            end
            for ia=[1 2 5 6],
                axis(o.ax(ia),'equal','tight');
            end
        end

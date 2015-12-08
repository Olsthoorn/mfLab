        function show(o,gr,iw,z,setlim,STCONC,debugLevel)
        % show(o,gr,iw,z[,STCONC|[] [,setlim]])
            
            if nargin < 4
                error('%s: use animate3D_TS_Obj.show(iwell,z)',mfilename);
            end
            
            % use Delta T and Delta C if nargin ==5;
            if  ~isempty(STCONC);
                o.STTEMP = STCONC{1};
                o.STCONC = STCONC{2};
                
                T_ = o.DT;
                C_ = o.DC;
                trange_ = o.dtrange;
                crange_ = o.dcrange;
            else
                C_ = o.C;
                T_ = o.T;
                trange_ = o.trange;
                crange_ = o.crange;
            end
            
            if exist('debugLevel','var')
                figname = sprintf('%s.Show, debugLevel = %d',mfilename,debugLevel);
            else
                figname = sprintf('%s.Show',mfilename);
            end
            figure('Name',figname,'position',get(0,'screensize'));
            grey = get(gcf,'color');
            
            o = o.setAxes(gr,setlim);

            time = o.well(1).t(2:end);
            

            %% Set up movie
            vidObj = VideoWriter([o.basename '.avi']);
            vidObj.FrameRate = o.framerate;
            vidObj.Quality   = o.quality;
            vidObj.open();

            %ix = round(median([o.well.ix]));
            try
                iy = hit(gr.yGr,o.well(iw).y);
                iz1 = hit(gr.zGr,z(  1));
                iz2 = hit(gr.zGr,z(end));
            catch %#ok
                error('%s: use animate3D_TS_Obj.show(iwell,z)',mfilename);
            end

            %% make movie
            for it=1:length(time)
                tts1 = [o.titleStr 'temp, t = ',datestr(time(it),'yyyy/mmm/dd')];
                tts2 = [o.titleStr 'conc, t = ',datestr(time(it),'yyyy/mmm/dd')];
                if it==1
                    ht1 = title(o.ax(1),tts1); set(o.ax(1),'clim',trange_([1 end]));
                    ht2 = title(o.ax(2),tts2); set(o.ax(2),'clim',crange_([1 end]));
                    ht3 = title(o.ax(3),tts1); set(o.ax(3),'clim',trange_([1 end]));
                    ht4 = title(o.ax(4),tts2); set(o.ax(4),'clim',crange_([1 end]));
                    ht5 = title(o.ax(5),tts1); set(o.ax(5),'clim',trange_([1 end]));
                    ht6 = title(o.ax(6),tts2); set(o.ax(6),'clim',crange_([1 end]));

                    
                    [~,hT2] = contourf(o.ax(3),gr.xc,gr.zc,XS(T_(it).values(iy,:,:  )),trange_,'edgecolor','none'); colorbar;
                    [~,hC2] = contourf(o.ax(4),gr.xc,gr.zc,XS(C_(it).values(iy,:,:  )),crange_,'edgecolor','none'); colorbar;
                    
                    % stream lines
                    if gr.Ny==1,
                        [~,hPT] = contour(o.ax(3),gr.xp,gr.zp,o.B(it).Psi,o.prange,'color',grey);
                        [~,hPC] = contour(o.ax(4),gr.xp,gr.zp,o.B(it).Psi,o.prange,'color',grey);
                    end
                    
                    if gr.Ny>1
                        [~,hT1] = contourf(o.ax(1),gr.xc,gr.yc,   T_(it).values(:,: ,iz1) ,trange_,'edgecolor','none'); colorbar;
                        [~,hC1] = contourf(o.ax(2),gr.xc,gr.yc,   C_(it).values(:,: ,iz1) ,crange_,'edgecolor','none'); colorbar;
                        [~,hT3] = contourf(o.ax(5),gr.xc,gr.yc,   T_(it).values(:,: ,iz2) ,trange_,'edgecolor','none'); colorbar;
                        [~,hC3] = contourf(o.ax(6),gr.xc,gr.yc,   C_(it).values(:,: ,iz2) ,crange_,'edgecolor','none'); colorbar;

%                         for ia=[1 2 5 6],
%                             axis(o.ax(ia),'equal','tight');
%                             set(o.ax(ia),'xlim',o.xLim,'ylim',o.yLim);
%                         end

                        well1 = well.plot3D(o.ax(1));
                        well2 = well.plot3D(o.ax(2));
                        well5 = well.plot3D(o.ax(5));
                        well6 = well.plot3D(o.ax(6)); % different wells on different axes

                    end
                    
                    set([o.ax(3) o.ax(4)],'ylim',[min(min(gr.ZBlay(:,:,end))),max(max(gr.ZTlay(:,:,1)))]);

                    % different wells because of different axes with their own handles

                    well3 = gr.plotXS(o.ax(3),[o.well.iy],o.well);
                    well4 = gr.plotXS(o.ax(4),[o.well.iy],o.well); % different wells on different axes

                else
                    set(ht1,'string',tts1);
                    set(ht2,'string',tts2);
                    set(ht3,'string',tts1);
                    set(ht4,'string',tts2);
                    set(ht5,'string',tts1);
                    set(ht6,'string',tts2);
                    
                    if gr.Ny>1
                        set(hT1,'zdata' ,   T_(it).values(: ,:,iz1));
                        set(hC1,'zdata' ,   C_(it).values(: ,:,iz1));
                        set(hT3,'zdata' ,   T_(it).values(: ,:,iz2));
                        set(hC3,'zdata' ,   C_(it).values(: ,:,iz2));
                        well1 = well1.plot3D(it);
                        well2 = well2.plot3D(it);
                        well5 = well5.plot3D(it);
                        well6 = well6.plot3D(it);
                    end

                    set(hT2,'zdata' ,XS(T_(it).values(iy,:,:  )));
                    set(hC2,'zdata' ,XS(C_(it).values(iy,:,:  )));

                    if gr.Ny==1
                        set(hPT,'zdata',o.B(it).Psi);
                        set(hPC,'zdata',o.B(it).Psi);
                    end

                    well3 = gr.plotXS(it,well3);
                    well4 = gr.plotXS(it,well4);
                end

                vidObj.writeVideo(getframe(gcf));

            end

            vidObj.close();

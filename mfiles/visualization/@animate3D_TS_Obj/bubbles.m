function bubbles(o,gr)
% animate3D_TS_Obj.bubbles()
%% Show DT (delta temp) and DC (delta conc) as bubbles
% TO 120908

% Get average injection temperature of the wells:
tM=Inf; tM=-Inf;
for iw=1:length(o.well)
    tm = min(tm,o.well(iw).UserData.tinj);
    tM = max(tM.o.well(iw).UserData.tinj);
end

tmean = 0.5*(tm+tM);
%% Axes positions

% ax1 = subplot(2,1,1,'xgrid','on','ygrid','on','zgrid','on',...
%         'nextplot','add','clim',o.dtrange([1 end]),'xlim',xlim,'ylim',ylim);
ax1 = axis('xgrid','on','ygrid','on','zgrid','on',...
        'nextplot','add','clim',o.dtrange([1 end]),'xlim',xlim,'ylim',ylim);
xlabel(ax1,'x [m]'); ylabel(ax1,'y [m]'); zlabel(ax1,'z [m]');
view(3);

% ax2 = subplot(2,1,2,'xgrid','on','ygrid','on','zgrid','on',...
%         'nextplot','add','clim',o.dcrange([1 end]),'xlim',xlim,'ylim',ylim);
% xlabel(ax2,'x [m]'); ylabel(ax2,'y [m]'); zlabel(ax2,'z [m]');
% view(3);

%% Addition by Jos Beemster 120820
% hlink = linkprop([ax1 ax2],{'CameraPosition','CameraUpVector'});
key = 'graphics_linkprop';
setappdata(ax1,key,hlink); 

time = o.well(1).t(2:end);

%% Set up movie
vidObj = VideoWriter([o.basename '.avi']);
vidObj.FrameRate = 3;
vidObj.Quality   = 80;
vidObj.open();

%% Make movie
for it=1:length(time)
    tts1 = [o.titleStr,', Temp [ C ], isosurface ',num2str(o.dtrange),', t=',datestr(time(it),'mmm-yy')];
%     tts2 = [o.titleStr,', Conc [g/L], isosurface ',num2str(o.dcrange),', t=',datestr(time(it),'mmm-yy')];
    if it==1
        ht1 = title(ax1,tts1);
%         ht2 = title(ax2,tts2);
        
%         axis ([ax1 ax2], 'equal');
%         axis ([ax1 ax2], 'auto z');
             
        well1 = gr.plot3D(ax1,o.well,'color','b');
%         well2 = gr.plot3D(ax2,o.well,'color','b');
    else
        delete(pt);
%         delete(pc);
        
        set(ht1,'string',tts1);
%         set(ht2,'string',tts2);
        
        well1 = gr.plot3D(it,well1);
%         well2 = gr.plot3D(it,well2);
    end
    
    for idt = length(o.dtrange):-1:1
        if o.dtrange(idt)<tmean,color='b'; else color='r'; end
        pt(idt) = patch(isosurface(gr.XM,gr.YM,gr.ZM,o.DT(it).values,o.dtrange(idt)),'parent',ax1);
        isonormals(gr.XM,gr.YM,gr.ZM,o.DT(it).values,pt(idt));
        set(pt(idt),'facealpha',0.2,'facecolor',color,'edgecolor','none');
        camlight; lighting('gouraud');
    end
    
%     for idc= length(o.dcrange):-1:1
%         if o.dcrange(idt)<cmean,color='b'; else color='r'; end
%         pc(idc) = patch(isosurface(gr.XM,gr.YM,gr.ZM,o.DC(it).values,o.dcrange(idc)),'parent',ax2);
%         isonormals(gr.XM,gr.YM,gr.ZM,o.DC(it).values,pc(idc));
%         set(pc(idc),'facealpha',0.2,'facecolor',color,'edgecolor','none');
%         camlight; lighting('gouraud');
%     end
        
        
    vidObj.writeVideo(getframe(gcf));

end

vidObj.close();

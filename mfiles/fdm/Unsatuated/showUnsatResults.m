
zm=0.5*(z(1:end-1)+z(2:end)); dz=abs(diff(z));

bodemCrop = ['bodem ',sprintf('%s,',bodem{:}), ' crop ' sprintf('%s,',crop{:})];

%% Show results
defaults = {'nextPlot','add','xGrid','on','yGrid','on','fontSize',12};
ylim = z([end,1]);

if showOut==true
    figure('position',[37   319   560   420]);
    ax1 = axes(defaults{:},'ylim',ylim,'xlim',[-35, 5],'xScale','lin');
    ttl1 = sprintf('h profile, %s time = %%5.0f d',bodemCrop);
    xlabel(ax1,'h [m]'); ylabel(ax1,'z [m]');

    figure('position',[690   316   560   420]);
    ax2 = axes(defaults{:},'ylim',ylim,'xlim',[0 1],'xScale','lin');
    ttl2 = sprintf('theta profile, %s time = %%5.0f d',bodemCrop);
    xlabel(ax2,'theta [m]'); ylabel(ax2,'z [m]');
    
    figure('position',[690    14   560   420]);
    ax3 = axes(defaults{:},'ylim',ylim,'xlim',[0 1],'xScale','lin');
    ttl3 = sprintf('aCrop profile, %s time = %%5.0f d',bodemCrop);
    xlabel('h [m]'); ylabel('z [m]');
    
    figure('position',[57    16   560   420]);
    ax4 = axes(defaults{:},'ylim',ylim,...
        'xlim',[tstart tend]);
    ttl4 = sprintf('head in Cover layer, %s time = %%5.0f d',bodemCrop);
    xlabel('time'); ylabel('head [m NAP]');
    
    set(ax4,'color','none');

    plot(ax1,h(:,1),zm(:),'r');      % initial pressure head
    plot(ax2,Theta(:,1),zm(:),'r');  % initial theta
    plot(ax3,Acrop(:,1),zm(:),'r');  % initial crop ET reduction factor
    datetick(ax4);
 
    for it=1:size(TPE,1)
        tttl1 = sprintf(ttl1,TPE(it,1)-TPE(1,1)+1); % days since start
        if it==1
            ht1 = title(ax1,tttl1);
            hp1 = plot(ax1,h(:,it),zm(:));
        else
            set(ht1,'string',tttl1);
            set(hp1,'xData',h(:,it));
        end

        tttl2 = sprintf(ttl2,TPE(it,1)-TPE(1,1)+1); % days since start
        if it==1
            ht2 = title(ax2,tttl2);
            hp2 = plot(ax2,Theta(:,it),zm(:));
        else
            set(ht2,'string',tttl2);
            set(hp2,'xData',Theta(:,it));
        end
        
        tttl3 = sprintf(ttl3,TPE(it,1)-TPE(1,1)+1); % days since start
        if it==1
            ht3 = title(ax3,tttl3);
            hp3 =plot(ax3,Acrop(:,it),zm(:));
        else
            set(ht3,'string',tttl3);
            set(hp3,'xData',Acrop(:,it));
        end
        
        tttl4 = sprintf(ttl4,TPE(it,1)-TPE(1,1)+1); % days since start
        if it==1
            ht4 = title(ax4,tttl4);
            hp4a = plot(ax4,TPE(it,1),h(end,it)+zm(end),'b','lineWidth',2);
            hp4b = plot(ax4,TPE(it,1),h(  1,it)+zm(end),'r','lineWidth',2);
        else
            set(ht4,'string',tttl4);
            set(hp4a,'xData',TPE(1:it,1),'yData',h(end,1:it)+zm(end));
            set(hp4b,'xData',TPE(1:it,1),'yData',h(  1,1:it)+zm(  1));
        end
        pause(0.05);
        
    end
    legend(ax1,'h(t1)','h(tEnd)');
    xlim = get(ax1,'xlim');
    text(xlim(1)+0.25*diff(xlim),mean(z),...
        sprintf('h(end)-h(1) = %.3f m',mean(h(:,end)-h(:,1))),'parent',ax1);
    
    legend(ax2,'theta(t1)','thetta(tEnd)');
    xlim = get(ax2,'xlim');
    text(xlim(1)+0.25*diff(xlim),mean(z),...
        sprintf('theta(end)-theta(1) = %.3f m',sum(dz.*(Theta(:,end)-Theta(:,1)))),'parent',ax2);

    legend(ax3,'acrop(t1)','acrop(tEnd)');
    xlim = get(ax3,'xlim');
    text(xlim(1)+0.25*diff(xlim),mean(z),'acrop');

    %drawnow;
end


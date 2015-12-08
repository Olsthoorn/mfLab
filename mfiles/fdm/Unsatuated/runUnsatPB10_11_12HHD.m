% Run unsat1
%% Determine what to show
showMeteo = true;
showOut   = true;
showBudget= true;
uitzakking= true;

%% get Meteo
meteoDir = '/Users/Theo/GRWMODELS/mflab/examples/TimeSeriesAnalysis/KNMI-data-files/';
load([meteoDir 'TPE_Rotterdam344']);  % loads TPS for Rotterdam

%% select simulation period
startDate = datenum(2009,10,1);
endDate   = datenum(2010,12,31);

TPE = TPE(TPE(:,1)>=startDate & TPE(:,1)<endDate,:);

if uitzakking==true
    tstart = datenum(2010,5,5);
    tend   = datenum(2010,8,1);
    I= TPE(:,1)>=tstart  & TPE(:,1)<= tend;
    TPE(I,2)=mean(TPE(I,2));
    TPE(I,3)=mean(TPE(I,3));
    str = sprintf('Van %s tot %s, %.0f d: mean P E P-E = %.4f %.4f %.4f, sum P E P-E = %.4f %.4f %.4f\n',...
        datestr(tstart),datestr(tend),tend-tstart,...
        mean(TPE(I,2)),mean(TPE(I,3)),mean(TPE(I,2)-TPE(I,3)),...
        sum( TPE(I,2)),sum( TPE(I,3)),sum( TPE(I,2)-TPE(I,3)));
    fprintf(str);
end

%TPE(:,3) = 0;

%% meteo Rotterdam
if showMeteo==true
    figure; axes('nextPlot','add','xGrid','on','yGrid','on','fontSize',12);
    xlabel('2010'); ylabel('m/d'); title('P and E in Rotterdam, station 344');

    
    bar(TPE(:,1),TPE(:,2),'b');
    plot(TPE(:,1),TPE(:,3),'r');
    datetick;
    
    % cumulative meteo Rotterdam vanaf 1 april 2010
    I = TPE(:,1)>=datenum(2010,4,1);
    figure; plot(TPE(I,1),cumsum(TPE(I,2)),'b',...
            TPE(I,1),cumsum(TPE(I,3)),'r',...
            TPE(I,1),cumsum(TPE(I,3)-TPE(I,2)),'k');
        datetick;
    xlabel('2010'); ylabel('m'); title('P and E en E-P cumulatief in Rotterdam, station 344');
    legend('P','E','E-P');
    grid on;
end

%% Run the model

[h,Qfix,Qnod,Qsto,QCau,Qz,z,bodem,crop,SS,Theta,Kh,Acrop]=unsat1(TPE,'ClapHornberger78SoilData');

zm=0.5*(z(1:end-1)+z(2:end)); dz=abs(diff(z));

bodemCrop = ['bodem ',sprintf('%s,',bodem{:}), ' crop ' sprintf('%s,',crop{:})];

%% Show results
defaults = {'nextPlot','add','xGrid','on','yGrid','on','fontSize',12};

if showOut==true
    figure('position',[37   319   560   420]);
    ax1 = axes(defaults{:},'ylim',z([end,1]),'xlim',[-5, 5],'xScale','lin');
    ttl1 = sprintf('h profile, %s time = %%5.0f d',bodemCrop);
    xlabel(ax1,'h [m]'); ylabel(ax1,'z [m]');

    figure('position',[690   316   560   420]);
    ax2 = axes(defaults{:},'ylim',z([end,1]),'xlim',[0 1],'xScale','lin');
    ttl2 = sprintf('theta profile, %s time = %%5.0f d',bodemCrop);
    xlabel(ax2,'theta [m]'); ylabel(ax2,'z [m]');
    
    figure('position',[690    14   560   420]);
    ax3 = axes(defaults{:},'ylim',z([end,1]),'xlim',[0 1],'xScale','lin');
    ttl3 = sprintf('aCrop profile, %s time = %%5.0f d',bodemCrop);
    xlabel('h [m]'); ylabel('z [m]');
    
    figure('position',[57    16   560   420]);
    ax4 = axes(defaults{:},'ylim',[-6 -3.5],...
        'xlim',[datenum(2010,1,1) endDate]);
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
            hp4 = plot(ax4,TPE(it,1),h(end,it)+zm(end),'lineWidth',2);
        else
            set(ht4,'string',tttl4);
            set(hp4,'xData',TPE(1:it,1),'yData',h(end,1:it)+zm(end));
        end
        pause(0.01);
        
    end
    legend(ax1,'h(end)','h(1)');
    xlim = get(ax1,'xlim');
    text(xlim(1)+0.25*diff(xlim),mean(z),...
        sprintf('h(end)-h(1) = %.3f m',mean(h(:,end)-h(:,1))),'parent',ax1);
    
    legend(ax2,'theta(end)','thetta(1)');
    xlim = get(ax2,'xlim');
    text(xlim(1)+0.25*diff(xlim),mean(z),...
        sprintf('theta(end)-theta(1) = %.3f m',sum(dz.*(Theta(:,end)-Theta(:,1)))),'parent',ax2);

    legend(ax3,'acrop(end)','acrop(1)');
    xlim = get(ax3,'xlim');
    text(xlim(1)+0.25*diff(xlim),mean(z),'acrop');

    %drawnow;
end

%% Check water balance
% Alles uit de berging
% Verify the storage coefficient
fprintf(' %11.4g',sum(SS .* (dz*ones(1,size(SS,2))))); fprintf \n
fprintf(' %11.4g',sum(diff(Theta,1,2)./diff(h,1,2)));  fprintf \n

if showBudget == true
    unsatBudget;

    fprintf \n\n\n
    fprintf('Water budget %s model unsat1\n', bodemCrop);
    fprintf('====================================================\n');
    fprintf('Total sum over all each of the individual water budget terms\n');
    fprintf(' %11s','sum(P-E)','sum(Qfix(:))','sum(Qsto(:))','sum(QCau(:))','sum(Qnod(:))','sum((theta(:,1)-theta(:,end)).*dz)');
    fprintf \n
    fprintf(' %11.4g',sum(TPE(:,2)-TPE(:,3)),sum(Qfix(:)),sum(Qsto(:)),sum(QCau(:)),sum(Qnod(:)),sum((Theta(:,1)-Theta(:,end)).*dz));
    fprintf \n\n

    fprintf(['Qsto must equal the difference in water content between start and finish\n',...
             'at least approximately, because of non-linearity:\n']);
    fprintf(' %11s','sum(Qsto(:))','sum((theta(:,end)-theta(:,1)).*dz');
    fprintf \n
    fprintf(' %11.4g',sum(Qsto(:)),sum((Theta(:,end)-Theta(:,1)).*dz));
    fprintf \n\n
    fprintf('The fixed and boundary inflows must equal the amount stored:\n');
    fprintf(' %11s','sum(Qfix(:)+QCau(:))','sum(Qsto(:))');
    fprintf \n
    fprintf(' %11.4g',sum(Qfix(:)+QCau(:)),sum(Qsto(:)));
    fprintf \n\n
end

if exist('str','var'), fprintf(str); end
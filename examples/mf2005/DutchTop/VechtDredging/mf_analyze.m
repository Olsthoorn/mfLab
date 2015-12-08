%Vechtdredging results
% TO 141209, use grid whereever possible (MODFLOW)

load underneath

H = readDat([basename '.HDS']);
B = readBud([basename '.BGT']);
B = mf_Psi(gr,B);

%% plot values over the full width of the cross section

ttl  = {'Prescribed head and head at NAP = -10 m'; ...
          'Total horizontal discharge in aquifer, eastward positive'; ...
          'Seepage through cover layer (upward positive)'};
leg  = {'Ground surface','Prescribed head','h2','h3','etc'};

%% Total cross section
    
ax1 = Conf.plotOverview(gr,H,B,ttl,leg);

%% Detail: same for cross section of viewWidth/2 m beyond dikes along river
viewWidth = 100;
X1=Conf.xL(strmatchi('VechtW',Conf.zoneNames))-viewWidth;
X2=Conf.xR(strmatchi('VechtE',Conf.zoneNames))+viewWidth;

ax2 = Conf.plotOverview(gr,H,B,ttl,leg,'xlim',[X1 X2]);

%% Plot stream-line patterns with heads and flows
ax3 = Conf.plotContours(gr,H,B,'lineWidth',2,'color','r');

%% Zoom in
set(ax3,'xlim',[X1 X2],'ylim',[-10 2]);

%% Compute and plot ground-bursting safety in all cells
% * total pressure (grains + water)
% * water pressure

%% Total pressure at the bottom of all cells
%
% $$ p_z = \int _{top} ^{bot} {\rho g dz} $$
%

g     =9.81;     % gravity      [N/kg]
rhow  =1000;     % density      [kg/m3]

WET = XS(WET);

p = cumsum (g * XS(gr.DZ).* (Conf.array2D('rhodry',gr.xGr) .*(1-WET) +  ...
                             Conf.array2D('rhowet',gr.xGr) .* WET) , 1);

% Water pressure at ground surface (model top).
% This is relevant in river where bottom is ground surface.
p0 = g * rhow * (H.values(1,:,1)-gr.Z(1,:,1)); p0(p0<0)=0;
p0 = ones(gr.Nz,1)*XS(p0);

% Total pressure including effect of water above ground surface
p=p+p0;

% Vertical specific discharge at bottom of each cell is taken for the
% specfic vertical discharge inside the bottom half of each cell. We need
% this because this changes the total pressure at the cell bottem relative
% to that in the cell center. Remember we compute the safety at cell
% bottoms, because its at cell bottoms where material properties suddenly
% change and bursting of soil may occur in case of upwelling with sand
% under clay.
Qz = B.term{strmatchi('FLOWL',B(end).label)};

%% Water pressure relative to cell bottom or zero if head below cell bottom 
sigmaW = rhow * g .* XS(( H(end).values-gr.ZBlay - (Qz./gr.DX).*(gr.DZ/2)./VK ));
sigmaW(sigmaW<0)=0;

%% Safety factor equals p/sigmaw, total pressure over water pressure at cell bottoms

% Safety factor = total pressure/water pressure
safety=p./sigmaW;

% Handle exeptions:
% 1) Water is always safe
safety(Conf.array2D('matindex',gr.xGr)==strmatchi('Water',Conf.matNames,'exact'))=Inf;

% 2) Dry cells are always safe
safety(~XS(WET)) = Inf;  % head below cell bottom

% 3) Downward flow always safe
safety(qy<0) = Inf;

% Cell bottom center coordinates, which is where the safety was computed
xcm = XS(gr.XM(1,:,:));
zcm = XS(gr.ZM(1,:,:));

%% Add safety status to plot, by plotting colored markers at all cell bottoms if any risk exists
ls=safety>=1.1 & safety<1.2;  plot(ax3(2),xcm(ls),zcm(ls),'b.','markersize',7); % small risk
ls=safety>=1.0 & safety<1.1;  plot(ax3(2),xcm(ls),zcm(ls),'m.','markersize',7); % medium risk
ls=safety< 1.0;               plot(ax3(2),xcm(ls),zcm(ls),'r.','markersize',7); % large risk

legend('head','stream lines','ground surface','prescribed head','small risk','medium risk','large risk');

%% Compute and show minimum safety factor in verticals of cross section
% Draw an extra plot showing the lowest safety in any vertical column. This
% is a graph of the computed safetyfactor versus x. This graph will be
% colored to indicate safe or risk at any x of the cross section.

% minimum safety along vertical lines
safetymin=min(safety,[],1);

% plot each range of safety in its own color
figure;
args = {'nextplot','add','xgrid','on','ygrid','on'};

ax4(2)=subplot(2,1,2,args{:});

xlabel(ax4(2),'x [m]');  ylabel(ax4(2),'safety factor [-]');
title(ax4(2),'Minimum safety factor against vertical soil breakup');

s=safetymin; ls=s < 1.0;         s(~ls)=NaN;  plot(ax4(2),gr.xm,s,'r','linewidth',2.0);
s=safetymin; ls=s >=1.0 & s<1.1; s(~ls)=NaN;  plot(ax4(2),gr.xm,s,'m','linewidth',1.5);
s=safetymin; ls=s >=1.1 & s<1.2; s(~ls)=NaN;  plot(ax4(2),gr.xm,s,'b','linewidth',1.0);
s=safetymin; ls=s >=1.2;         s(~ls)=NaN;  plot(ax4(2),gr.xm,s,'g','linewidth',0.7);

set(ax4(2),'xlim',xlim);

legend(ax4(2),'unsafe (<1.0)','large risc (1.0-1.1)','risc (1.1-1.2)','safe (>1.2)');

% Couple the x-axis with those of the previous figure (which has to axes
% ax7a (material color patches) and ax7 the contours and lines.
hlink1 = linkprop([ax3,ax4(2)],'xlim');

% To show Detail of dredging set detail=1 else set detail=0.
% Overrule existing xlim (and ylim) to show more detail
detail = 1;
if detail
    set(ax4(2),'xlim',0.5*viewWidth*[-1 1]);
else
    set(ax4(2),'xlim',gr.xGr([1 end]),'ylim',5*[floor(min(Z(:)/5)) ceil(max(Z(:)/5))]); %#ok<*UNRCH>
end

ax4(1) = subplot(2,1,1,args{:});
xlabel('x [m]'); ylabel('elevation [m]'); title('Cross section ground surface and heads');

plot(ax4(1),gr.xm,Conf.array2D('top'    ,xGr),'g','linewidth',2);
plot(ax4(1),gr.xm,XS(H(end).values(1,:,1:4)));
legend('Ground Surface','prescribed head','head(2)','head(3)','head(4)');

hlink2 = linkprop(ax4,'xlim');
%% Water budget
fprintf('Total infiltration through river bottom is %.2f m2/d\n',...
    sum(Qz(gr.xm>=X1 & gr.xm<=X2).*gr.dx(gr.xm>=X1 & gr.xm<=X2)));

%% Overview of config and material properties
Conf.show('config');
Conf.show('materials');

%% Check safety computation (verification)
% This code prints a table for the center of  one of the zones depending on
% iCase and cases below. This table  can be used to verify the model by
% doing hand calculations easily.

cases = {'VechtW','VechtE','Polder1W','PolderW'};
% 
iCase =0;

if iCase
    i = round(median(find(gr.xm>Conf.xL(strmatchi(cases{iCase},Conf.zoneNames)) &...
                          gr.xm<Conf.xR(strmatchi(cases{iCase},Conf.zoneNames)))));
    % print header of table
    fprintf(' j');
    fprintf('%13s','Z','DZ','WET','Phi','Qy','dPhiQ','HK','VK','p','sigmaW','safety');
    fprintf('\n');

    % specific vertical discharge at selected column for all cell bottoms
    qv = [Qy(:,i);0]/Dx(i);

    % print table
    for j=1:Nz
        fprintf(['%2d ' repmat(' %12g',[1,11]) '\n'],...
            j,Z(j+1,i),DZ(j,i),WET(j,i),Phi(j,i),qv(j),qv(j)*DZ(j,i)/2/VK(j,i),HK(j,i),VK(j,i),p(j,i),sigmaW(j,i),safety(j,i));
    end
end                   
                   


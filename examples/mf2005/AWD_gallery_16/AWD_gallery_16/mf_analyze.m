%% ==== Visualize output ========

verbose=1;

load name; load(basename);
load underneath;

%% Reading the data
ixs=[30,grid.Nx-30];

H=maskHC(readDat([basename, '.hds']),1000);
B=       readBud([basename, '.bgt']);
B=mf_Psi(B,ixs,'y');

%%

plotmesh;
    
for it=1:length(B) 
    t=H(it).totim;
    if it==1
        [~,hdl]=contour(grid.xm,grid.ym,H(it).values(:,:,iLay),hrange); 
        htl=title(sprintf('t=%d d, Layer %d, z=[%.2f %.2f] mNAP, dPhi=%.2f m',t,iLay,zGr([iLay iLay+1]),dphi));
    else
        set(htl,'string',sprintf('Layer %d, z=[%.2f %.2f] mNAP, dPhi=%.2f m',t,iLay,zGr([iLay iLay+1]),dphi))
        set(hdl,'zdata',H(it).values(:,:,iLay)); 
    end
    drawnow;
end

%% Dealing with the budget

QCH=B(end).term{strmatchi('CONSTANTHEAD',B(end).label)};

zonebudget(B);

zones=unique(abs(IBOUND(:)));
for i=1:length(zones)
    zonebudget(B,abs(IBOUND),zones(i));
end

%% plot Q in Drain

% Collect the flow along each drain in a separate array
figure
Q=[]; NDr=length(Drain); Nt=length(B);
for i=NDr:-1:1
    DX{i}=grid.DX(Drain(i).I);
    Q(i).A=NaN(Nt,numel(Drain(i).I));
    for it=1:Nt
        Q (i).A(it,:)=B(it).term{strmatchi('CONSTANTH',B(it).label)}(Drain(i).I);
    end
    QS(i).A=cumsum(Q(i).A,2);
    q( i).A=Q(i).A./(ones(Nt,1)*DX{i});
    subplot(2,1,i);
    plot(Drain(i).L,QS(i).A);
    title(sprintf('Cumulative extraction by drain from toe to heel (x=0) in {m}^{3}/d, drain %d',i));
end

%% Comulative flow along the drain

% first mark the piezometers that are either in the Drain P, or just above A B or C

CPIPES={'C105','C106','C107','C108','C109','C110','C111','C112','C113','C114','C115','C116'};

CP={'C105','C106','C107','C108','C109','C110'};
CA={'C105M_A','C106M_A','C107M_A','C108M_A','C109M_A','C110M_A'};
CB={'C105M_B','C106M_B','C107M_B','C108M_B','C109M_B','C110M_B'};
CC={'C105M_C','C106M_C','C107M_C','C108M_C','C109M_C','C110M_C'};

for i=1:length(Piez), Piez(i).type='none'; end
for i=1:length(CP),
    I=strmatchi(CP{i},{Piez.name},'exact');
    if ~isempty(I)
        Piez(I).type='p';
    end
end
for i=1:length(CA),
    I=strmatchi(CA{i},{Piez.name},'exact');
    if ~isempty(I),
        Piez(I).type='a';
    end
end
for i=1:length(CB),
    I=strmatchi(CB{i},{Piez.name},'exact');
    if ~isempty(I),
        Piez(I).type='b';
    end
end
for i=1:length(CC),
    I=strmatchi(CC{i},{Piez.name},'exact');
    if ~isempty(I),
        Piez(I).type='c';
    end
end

%% second: add the computed head to the Piez struct

for i=1:length(Piez)
    [Piez(i).ix,Piez(i).iy,Piez(i).iz]=xyzIndex(Piez(i).x,Piez(i).y,Piez(i).z_bot,xGr,yGr,zGr);
    
    try
        Piez(i).I =cellIndex(Piez(i).ix,Piez(i).iy,Piez(i).iz,Nx,Ny,Nz);
        Piez(i).Phi=Phi(Piez(i).I);
    catch
       fprintf('Piez(%d), %s lies outside the grid\n',i,Piez(i).name);
    end
end

%% third: plot it

figure; hold on; grid on; xlabel('x [m]'); ylabel('y [m]'); grid on;

title(sprintf('DR16 %s; Q_{South}=%.0f, Q_{North}=%.0f Q_{tot}=%.0f m^3/d',...
    datestr(time,'dd-mmm-yyyy'),Drain(1).QT, Drain(2).QT, Drain(1).QT+Drain(2).QT));

% Get grid location of marked piezometers
ip=strmatchi('p',{Piez.type}); IP=[Piez(ip).I];  % first is in list of Piez, second global grid indices
ia=strmatchi('a',{Piez.type}); IA=[Piez(ia).I];
ib=strmatchi('b',{Piez.type}); IB=[Piez(ib).I];
ic=strmatchi('c',{Piez.type}); IC=[Piez(ic).I];

% Get the correct head for plotting
hp=NaN(size(ip)); 
for i=1:length(ip),
    h=Piez(ip(i)).head(Piez(ip(i)).t==time);
    if ~isnan(h) hp(i)=h; end
end
ha=NaN(size(ia));
for i=1:length(ia),
    h=Piez(ia(i)).head(Piez(ia(i)).t==time);
    if ~isnan(h), ha(i)=h; end
end
hb=NaN(size(ib));
for i=1:length(ib),
    h=Piez(ib(i)).head(Piez(ib(i)).t==time);
    if ~isnan(h), hb(i)=h; end
end
hc=NaN(size(ic));
for i=1:length(ic),
    h=Piez(ic(i)).head(Piez(ic(i)).t==time);
    if ~isnan(h), hc(i)=h; end
end


if 1
    plot(XM(IP),hp,'rp');  % inside the Drain
    plot(XM(IA),ha,'bp');  % just outside
    plot(XM(IB),hb,'gp');  % more outside
    plot(XM(IC),hc,'kp');  % still more outside (stagnation points)
else
    plot(XM(IP),Phi(IP),'rp');  % inside the Drain
    plot(XM(IA),Phi(IA),'bp');  % just outside
    plot(XM(IB),Phi(IB),'gp');  % more outside
    plot(XM(IC),Phi(IC),'kp');  % still more outside (stagnation points)
end


%%
for i=1:length(Drain)
    Ip=cellIndex(Drain(i).ix,Piez(ip(1)).iy,Drain(   i).iz,Nx,Ny,Nz);
    Ia=cellIndex(Drain(i).ix,Piez(ia(1)).iy,Piez(ia(1)).iz,Nx,Ny,Nz);
    Ib=cellIndex(Drain(i).ix,Piez(ib(1)).iy,Piez(ib(1)).iz,Nx,Ny,Nz);
    Ic=cellIndex(Drain(i).ix,Piez(ic(1)).iy,Piez(ic(1)).iz,Nx,Ny,Nz);

    plot(XM(Ip),Phi(Ip),'r');
    plot(XM(Ia),Phi(Ia),'b');
    plot(XM(Ib),Phi(Ib),'g');
    plot(XM(Ic),Phi(Ic),'k');

end

legend({'C-well','a (just outside)','b (more outside)','c (stagnation pnt)',...
    'Drain','computed at a','computed at b','computed at c'});

%%
% hdrange=contourRange(H,0.05);
% figure; hold on; xlabel('x [m]'); ylabel('elevation [NAP]'); title('cross section');
% contourf(ym,zm,squeeze(permute(H(end).values(:,Drain(1).ix(end),:),[3,1,2])),hdrange);
% set(gca,'xlim',[-16 16],'ylim',[zGr(end) zGr(1)]);
% 
% frf=permute(B(end).term{strmatchi('FLOWR',B(end).label)}(:,Drain(1).ix(end),:),[3,1,2]);
% Psi=flipud(cumsum(flipud([frf zeros(size(frf,1))]),1));
% 
%psirange=contourRange(Psi,0.5);

%contour(yGr(2:end-1),zGr,Psi(:,1:end-1),'y');

%% Following steps would be to start dealing ]with resistance in and outside
% of the Drain. Or to start using CDF but this would take several days to
% implement test and get reasonalble results with. There is no time to do


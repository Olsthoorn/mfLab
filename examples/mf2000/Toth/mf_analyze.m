%% Analyzing output of the model
% TO 091011 091129

load('name.mat') % get basename stored in file name.mat
load(basename);  % having retrieved baasename load the data in basename.mat
load underneath  % to get gr object

ayear = 365.24;
aweek = 7;

mf_checkdir;

H = readDat([basename,'','.hds']); % read the unformatted head file
B = readBud([basename '.BGT']);
B = mf_Psi(B);

hrange = ContourRange(H,50);
prange = ContourRange(B,50,'Psi');

defaults = {'nextplot','add','fontsize',12,'xlim',gr.xc([1 end]),'clim',hrange([1 end])};


figure('name','Example Toth','position',screenPos(0.75));

pos1 = [0.1 0.7 0.8 0.2];
pos2 = [0.1 0.1 0.8 0.5];

ax1 = axes('position',pos1,'ylim',[gr.zGr(1)-5 max(hdrn)+5],defaults{:},'xgrid','on','ygrid','on');
ax2 = axes('position',pos2,'ylim',gr.zc([end 1])           ,defaults{:});

xlabel(ax2,'x [m]');
ylabel(ax2,'elev [m]');
ylabel(ax1,'elev [m]');

plot(ax1,gr.xm,hdrn,'r');

%% Movie of cross section through time

tsa='Toth''s water table = %.1f weeks';
tsb='Toth''s flow system with heads and stream lines = %.1f weeks';

time = [H.totim]/aweek;

vidObj = VideoWriter(basename);
vidObj.FrameRate= 1;
vidObj.Quality = 80;
vidObj.open;

for it=1:length(time)
    
    ts1   =  sprintf(tsa,time(it));
    ts2   =  sprintf(tsb,time(it));
    
    if it==1        
         h = plot(ax1,gr.xm,hdrn,'r',gr.xm,H(1).values(1,:,1),'b');
        
          ht1=title(ax1,ts1);
          ht2=title(ax2,ts2);
        [~,hh] = contourf(ax2,gr.xc,gr.zc,XS(H(it).values),hrange,'edgecolor','none');
        [~,hp] = contour( gr.xp,gr.zp,B(it).Psi,prange,'color','b');
    else
        set( h,'ydata',H(it).values(1,:,1));
        set(ht1,'string',ts1);
        set(ht2,'string',ts2);
        set(hh,'zdata' ,XS(H(it).values));
        set(hp,'zData' ,B(it).Psi);
    end
    vidObj.writeVideo(getframe(gcf));
end

vidObj.close;

%% ZoneBudget

zonebudget(B)

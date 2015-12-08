% MF_ANALYZE: Visualize and interpret model output
%
% USAGE:
%   mf_analyze
%
% Analyzing output of the model
% TO 091011 091129 110429 120514

clear variables; close all;

load('name.mat') % loads name.mat, contains the variable "basename"                 
load(basename);  % having retrieved baasename load the data in basename.mat
load underneath;

ttl='A rainwater lense in action, streamlines vary continuously, but displacement is small ';

layer=3;
ayear=365.25;

%% Reading data
H=readDat([basename,'.hds']); % read the unformatted head file
C=readMT3D('MT3D001.UCN');    % read the unformatted head file
B=readBud([basename '.bgt']);
B=mf_Psi(B);

%% Contour Ranges
hrange = ContourRange(H,50);
crange = [0 ContourRange(C,50)];
prange = ContourRange(B,100,[],'Psi');

%% This is the figure that will receive all plot instructions
figure; 
grey=get(gcf,'color');

axpos = [0.1 0.1 0.8 0.8];

ax=axes('nextplot','add','position',axpos,'color',[0.6 0.6 1],...
    'xlim',gr.xc([1 end]),'ylim',gr.zc([end-1 1]),'clim',crange([1 end]));

ts1 = sprintf('Build-up of rainwater lens. Streamlines vary daily, but displacement is small. dPsi=%.2f m2/d, t= %%.1f y',min(diff(prange)));

%% Conc in section

vidObj=VideoWriter(basename);
vidObj.FrameRate=10;
vidObj.Quality=80;
vidObj.open();

time = [H.totim];

for it=1:length(time)
    
    ts2 = sprintf(ts1,time(it)/ayear);
    
    if it==1
        ht = title(ts2);
        hf = fill([gr.xc gr.xc(end:-1:1)],[gr.ZBlay(1,:,end) gr.ZTlay(1,end:-1:1,1)],'y');        
        [~,hc]=gr.contourf(ax,XS(C(it).values),crange,'EdgeColor','none');
        [~,hp]=gr.streamlines(ax,B(it).Psi,'color',grey);
        hpl = plot(ax,gr.xc,H(it).values(1,:,3),'b');
        
    else
        set(ht,'string',ts2);
        set(hc,'zdata',XS(C(it).values));        
        gr.streamlinesUpdate(hp ,B(it).Psi);
    end
    
    vidObj.writeVideo(getframe(gcf));

end

vidObj.close;

%% Zone budget
zonebudget(B);

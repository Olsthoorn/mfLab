function aviobj=movie3D(what,H,xGr,yGr,zGr,conts,color_switch_temp,well,tpvwjpg,tpvwxy,basename,film)
% aviobj=movie3D(what,H,xGr,yGr,zGr,conts,color_switch_temp,well,tpvwjpg,tpvwxy,film)
% what any text indicating what to plot is used in plot title
% make a 3D movie of the ATES simulation
% H conc struct (may also be head struct)
% x, y, z coordinates
% conts = isoplanes ot be drawn
% color_switch_temp above this temp patches are red below blue
% well well struct
% topvwjpg, a jpg file name with picture of topview
% tpvwxy coordinates of topview [xLL yLL xUR uUR]
% basename, basename for avi file
% film make avi if true
% TO 091121

[xGr,yGr,zGr,xm,ym,zm]=modelsize3(xGr,yGr,zGr);

if isempty(findstr('.jpg',tpvwjpg)), tpvwjpg=[tpvwjpg '.jpg']; end

topview = imread(tpvwjpg);
[topview,map]=rgb2ind(topview,256);
topview=topview(end:-1:1,:,:);

cols =H(1).cols;
rows =H(1).rows;
lays =H(1).lays;

xlim=sort(xGr(cols([1 end]))); if ~isempty(tpvwxy), xlim=tpvwxy([1 3]); end
ylim=sort(yGr(rows([1 end]))); if ~isempty(tpvwxy), ylim=tpvwxy([2 4]); end
zlim=sort(zGr(lays([1 end])));

conts=fliplr(sort(conts(:))');

figure
if film
    aviobj=avifile([basename ' 3dAVI'],...
        'compression','none',...
        'fps',8,...
        'colormap',colormap,...
        'videoname',...
        'KWO well');
end

for it=1:length(H)
    
    if isfield(H,'totim'), totim=H(it).totim; else totim=H(it).time; end

    title(sprintf('%s time=%.0f d',what,totim));
    
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');

    set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);
    set(gca,'XTicklabel',get(gca,'XTick'))
    set(gca,'YTicklabel',get(gca,'YTick'))

    if length(conts)>1, caxis(sort(conts([1 end]))); end
    
    grid('on'); hold('on');view(30,30);
    
    % plot the wells
    if it==1,
        for iW=1:length(well)
            if well(iW).isCold, clr='b'; else clr='r'; end
            plot3(well(iW).x([1,1]),well(iW).y([1,1]),[well(iW).scrtop well(iW).scrbot],...
                'color',clr,'linewidth',2);
            plot3(well(iW).x([1,1]),well(iW).y([1,1]),[       0        well(iW).scrtop],...
                'color','k','linewidth',1);
           % plot3(well(iW).x,well(iW).y,0,'marker','o','markerfacecolor','yellow','markeredgecolor','yellow')
        end
        
        if exist('topview','var')
            [m,n]=size(topview); dx=diff(tpvwxy([1 3]))/(n-1); dy=diff(tpvwxy([2 4]))/(m-1);

            [X,Y]=meshgrid(tpvwxy(1):dx:tpvwxy(3),tpvwxy(2):dy:tpvwxy(4));

            hsurf=surface(X,Y,zeros(m,n),topview,...
               'FaceColor','texturemap',...
               'EdgeColor','none',...
               'CDataMapping','direct',...
               'clipping','off',...
               'FaceAlpha',1,...
               'visible','on');
           colormap(map);
        end
    
        hp   =NaN(size(conts));
        hpold=NaN(size(conts));
        for ipatch=1:length(conts)
            if conts(ipatch)>color_switch_temp,fcolor='red'; else fcolor='blue'; end
            H(it).values(H(it).values<0)=NaN;
            hp(ipatch)=patch(isosurface(xm(cols),ym(rows),zm(lays),H(it).values,conts(ipatch)));
            isonormals(xm(cols),ym(rows),zm(lays),H(it).values,hp(ipatch));
            set(hp(ipatch),'FaceColor',fcolor,'EdgeColor','none',...
                 'FaceAlpha',0.5,...
                 'FaceVertexAlpha',0.5,...
                 'visible','off'); % hold on
             camlight('right');
             material('shiny');
             lighting('phong');
        end
    else
        for ipatch=1:length(conts)
            if conts(ipatch)>color_switch_temp,fcolor='red'; else fcolor='blue'; end
            H(it).values(H(it).values<0)=NaN;
            hpold(ipatch)=hp(ipatch);
            hp(ipatch)=patch(isosurface(xm(cols),ym(rows),zm(lays),H(it).values,conts(ipatch)));
            isonormals(xm(cols),ym(rows),zm(lays),H(it).values,hp(ipatch));
            set(hp(ipatch),'FaceColor',fcolor,'EdgeColor','none',...
                 'FaceAlpha',0.5,...
                 'FaceVertexAlpha',0.5,...
                 'visible','off'); % hold on
             delete(hpold(ipatch));
             set(hp(ipatch),'visible','on');
        end

    end

   if film, F=getframe(gcf); aviobj=addframe(aviobj,F); end
   
end

if film, aviobj=close(aviobj); end
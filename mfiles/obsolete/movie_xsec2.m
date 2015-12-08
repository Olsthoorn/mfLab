function aviobj=movie_xsec2(H,xGr,yGr,zGr,Sx,Sy,Sz,what,conts,color_switch_temp,jpegfile,coords,ZImage,basename,film)
% MOVIE_XSEC: Makes movie of vertical cross section of model using output H
%   where output may be any concentrationt struct of head struct read in
%   using readDAT and readMT3D
%
% USAGE:
%   movie_xsec(H,hpc,vpc,what,conts,dir,idx,basename,jpegfile,coords,film)
%   hpc is hor plot coordinate, vpc is vertical plot coordinate (either
%   x,y,z depending on dir
%   z  are cell center elevations for vertical plot
%  what is either HEADS or something else, the text will be used in labels
%   contour interval or list of contours
%   but to plot heads, what must be HEAD
%   dir is plortorientation: one of 'XZ', 'YZ', 'XY', 'YX'
%   idx is the row or color number for the contours
%   basename is used in title
%   film=boolean, if~=0 then make movie else just show figs
%
% TO 091121 091223

if isempty(H)
    error('Length data struct is zero, can''t visualize.\n');
end

[xGr,yGr,zGr,xm,ym,zm]=modelsize3(xGr,yGr,zGr);

[XM,YM,ZM]=meshgrid(xm(H(1).cols),ym(H(1).rows),zm(H(1).lays));


%% start making movie

if isempty(coords)
    xlim=xm([1,end]);   ylim=sort(ym([1 end])); zlim=sort(zm([1 end]));
else
    xlim=coords([1 3]); ylim=coords([2 4]);     zlim=sort(zm([1 end]));
end

dstr=sprintf('contours=%s',sprintf(' %g',conts));

figure; hold on; grid on;

set(gcf','doublebuffer','on');

set(gca,'xlimmode','manual',...
    'ylimmode','manual',...
    'zlimmode','manual',...
    'climmode','manual',...
    'alimmode','manual');


set(gca,'xlim',xlim,'ylim',ylim,'zlim',zlim);

if length(conts)>1, caxis(sort(conts([1 end]))); end

xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');

if film
    aviobj=avifile([basename ' vertical'],'compression','none',...
    'fps',8,'colormap',colormap,'videoname','KWO well');
end

clrmap2=colormap;  % get default colormap


for it=1:length(H)  % we simply plot all times in struct (no selection)    

    if it==1
        if ~isempty(jpegfile) && ~ isempty(coords)
            topview = imread(jpegfile); [topview,map]=rgb2ind(topview,256); topview=flipud(topview);

            [m,n]=size(topview); dx=diff(coords([1 3]))/(n-1); dy=diff(coords([2 4]))/(m-1);

            [X,Y]=meshgrid(coords(1):dx:coords(3),coords(2):dy:coords(4));

            colormap(map);
            surface(X,Y,ZImage*ones(m,n),topview,...
            'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
            brighten(0.8);
            clrmap1=colormap;
            
            colormap([clrmap1;clrmap2]); % combinede colormap image must be
            % first because it uses direct indexing 0-255, while the second
            % applies scaling of CData
            
            L1=length(clrmap1);
            L2=length(clrmap2);
            
            CLim=newclim(L1+1,L1+L2,conts(1),conts(end),L1+L2); % to correctly strmatchi in part of the colormap
            caxis(sort(CLim));
            %colorbar;  % see the colorbar to view both parts of colorbar

        end
        axis equal
        axis tight
        view(2);
        daspect([1,1,0.3])
    else
        hdlOld=hdl;
    end
    
    hdl=contourslice(XM,YM,ZM,H(it).values,Sx,Sy,Sz,conts);
    for i=1:length(hdl)
        ctemp=get(hdl(i),'FaceVertexCData');
        
        if ctemp(1)<color_switch_temp,
            set(hdl(i),'EdgeColor','b');
        else
            set(hdl(i),'EdgeColor','r');
        end
    end
    if exist('hdlOld','var')
        delete(hdlOld);
    end
    
    if isfield(H,'totim'), t=H(it).totim; else t=H(it).time ; end

    title(sprintf('%s time=%6.0f d; %s',what,t,dstr));

    if film,
        F=getframe(gcf); aviobj=addframe(aviobj,F); % close(gcf);
    end
end
if film, aviobj=close(aviobj); end


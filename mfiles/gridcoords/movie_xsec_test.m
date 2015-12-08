function aviobj=movie_xsec_test(H,gr,Sx,Sy,Sz,what,conts,color_switch_temp,jpegfile,coords,ZImage,basename)
%MOVIE_XSEC_TEST generate movie for vertical cross section of model using output H
%
% Example:
%     aviobj=movie_xsec(H,xGr,yGr,what,conts,color_switch_temp,jpegfile,coords,ZImage,basename,film)
% 
% where output may be any concentrationt struct of head struct read in
% using readDAT and readMT3D
% xGr,yGr,zGr, grid coordinates
% Sx,Sy,Sz, lines where the contourslices are required
% what is either HEADS or something else, the text will be used in labels
% contour interval or list of contours
% but to plot heads, what must be HEAD
% basename is used in title
% film=boolean, if~=0 then make movie else just show figs
%
% Used in:
%   mflab/examples/swt_v4/ATES-WKO/DamSqr-Ruben
%   mflab/examples/swt_v4/ATES-WKO/mfLabKWO-Dam/Dam_Pest/euronext
%   mflab/examples/swt_v4/ATES-WKO/mfLabKWO-Dam/Dam_Pest/euronext/metingen
%
% ToDo: replace by animateObj and modernize the rubish TO 130428)
%
% See also: animateObj
%
% TO 091121

if isempty(H)
    error('Length data struct is zero, can''t visualize.\n');
end

%% start making movie

if isempty(coords)
    xlim = gr.xm(  [1,end]);
    ylim = gr.ym(  [end 1]);
    zlim = gr.zLay([end 1]);
else
    xlim = coords([1 3]);
    ylim = coords([2 4]);
    zlim = gr.zm([1 end]);
end

dstr=sprintf('contours=%s',sprintf(' %g',conts));

figure('position',screenPos(0.75),'doublebuffer','on');

axes('nextplot','add','xGrid','on',yGrid','on','xlim',xlim,'ylim',ylim,'zlim',zlim);

xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');

if length(conts)>1, caxis(sort(conts([1 end]))); end

aviobj=avifile([basename ' vertical'],'compression','none',...
        'fps',8,'colormap',colormap,'videoname','KWO well');

clrmap2=colormap;  % get default colormap

for it=1:length(H)  % we simply plot all times in struct (no selection)    

    if it==1
        if ~isempty(jpegfile) && ~ isempty(coords)
            topview = imread(jpegfile);
            [m,n]=size(topview); dx=diff(coords([1 3]))/(n-1); dy=diff(coords([2 4]))/(m-1);
            topview=topview(end-1:-1:1,:,:);

            [X,Y]=meshgrid(coords(1):dx:coords(3),coords(2):dy:coords(4));
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
            %colorbar;  % if you line see the colorbar

        end
        axis equal
        axis tight
        view(2);
        daspect([1,1,0.3])
    else
        hdlOld=hdl;
    end
    
    hdl=contourslice(gr.XM,gr.YM,gr.ZM,H(it).values,Sx,Sy,Sz,conts);
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

    aviobj=addframe(aviobj,getframe(gcf));
end

aviobj=close(aviobj);


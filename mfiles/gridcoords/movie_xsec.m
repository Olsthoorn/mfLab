function movie_xsec(H,hpc,vpc,what,conts,dir,strmatchi,jpegfile,coords,basename,film)
%MOVIE_XSEC Makes movie of vertical cross section of model using output H
%   where output may be any concentrationt struct of head struct read in
%   using readDAT and readMT3D
%
% USAGE:
%   movie_xsec(H,hpc,vpc,what,conts,dir,idx,basename,jpegfile,coords,film)
%   hpc is hor plot coordinate, vpc is vertical plot coordinate (either
%   x,y,z depending on dir
%   z  are cell center elevations for vertical plot
%   what is either HEADS or something else, the text will be used in labels
%   contour interval or list of contours
%   but to plot heads, what must be HEAD
%   dir is plortorientation: one of 'XZ', 'YZ', 'XY', 'YX'
%   idx is the row or color number for the contours
%   basename is used in title
%   film=boolean, if~=0 then make movie else just show figs
%
% Used in:
%    mflab/examples/swt_v4/ATES-WKO/DamSqr-Ruben/Monowells
%    mflab/examples/swt_v4/ATES-WKO/mfLabKWO-Dam/Dam_Pest/euronext
%    mflab/examples/swt_v4/ATES-WKO/mfLabKWO-Dam/Dam_Pest/euronext/metingen
%
% ToDo: replace by animateObj and replace well by wellObj, modernize (TO 130428)
%
% SEE ALSO: movie3D
%
% TO 091121

if isempty(H)
    error('Length data struct is zero, can''t visualize.\n');
end

dir =upper(dir);

hpm = 0.5*(hpc(1:end-1)+hpc(2:end));  % coordinates along hor  axis
vpm = 0.5*(vpc(1:end-1)+vpc(2:end));  % coordinates along vert axis

%% determine hor coordinates depending on direction dir
switch dir
    case 'XZ'
        labelx='x [m]';
        labely='z [m]';
        xmplot=hpm(H(1).cols);
        ymplot=vpm(H(1).lays);            
    case 'YZ'
        labelx='y [m]';
        labely='z [m]';
        xmplot=hpm(H(1).rows);
        ymplot=vpm(H(1).lays);
    case 'XY'
        labelx='x [m]';
        labely='y [m]';
        xmplot=hpm(H(1).cols);
        ymplot=vpm(H(1).rows);
    case 'YX' 
        labelx='y [m]';
        labely='x [m]';
        xmplot=hpm(H(1).rows);
        ymplot=vpm(H(1).cols);
    otherwise
        error('Dir must be either ''XZ'' or ''YZ'' or ''XY'' or ''YX'', not <<%s>>',dir);
end

%% start making movie
figure; hold on;
grid on;

if isempty(coords)
    xlim=xmplot([1,end]); ylim=sort(ymplot([1 end]));
else
    xlim=coords([1 3]); ylim=coords([2 4]);
end

set(gca,'xlim',xlim,'ylim',ylim);
xlabel(labelx); ylabel(labely);

if film
    mov=avifile([basename ' vertical'],'compression','none',...
    'fps',8,'colormap',colormap,'videoname','KWO well');
end

if ~isempty(jpegfile) && ~ isempty(coords)
    topview = imread(jpegfile); [topview,map]=rgb2ind(topview,256); topview=flipud(topview);

    [m,n]=size(topview); dx=diff(coords([1 3]))/(n-1); dy=diff(coords([2 4]))/(m-1);

    [X,Y]=meshgrid(coords(1):dx:coords(3),coords(2):dy:coords(4));

    surface(X,Y,zeros(m,n),topview,...
    'FaceColor','texturemap',...
    'EdgeColor','none',...
    'CDataMapping','direct','clipping','off','FaceAlpha',1,'visible','on');
    colormap(map)
    
 %   alpha(0.3);
end

dstr=sprintf('contours=%s',sprintf(' %g',conts));

for i=1:length(H)  % we simply plot all times in struct (no selection)    
    if length(conts)>1, caxis(sort(conts([1 end]))); end
    h=turnselect(H(i).values,strmatchi,dir);
    if i==1
        [c,hdl]=contourf(xmplot,ymplot,h,conts);
        axis equal
        axis tight
        alpha(0.3);
 %      colorbar;
    else
        set(hdl,'Zdata',h);
    end

    if isfield(H,'trpstp'), sstp=sprintf('trpst=%d',H(i).trpstp);
    else                    sstp=sprintf('tstp=%d',H(i).tstp);
    end
    
    if isfield(H,'totim'), t=H(i).totim; end
    if isfield(H,'time') , t=H(i).time ; end

    title(sprintf('%s Stress period=%d; %s; time=%.3gd; %s',what,H(i).period,sstp,t,dstr));

    if film,
        F=getframe(gcf); mov=addframe(mov,F); % close(gcf);
    end
end
if film, mov=close(mov); end


function h=turnselect(values,idx,dir)
% h=turnselect(values,idx,dir)
% turn the matrix according to dir for plotting on the correct face idx
% return this 2D face in the correct position for plotting
% dir as above
% idx may be a vecot, in which case the average over idx is taken (no
% weighting of layer thickness nor porosities)

switch dir
        case 'XZ',
            h=permute(mean(values(idx,:,:),1),[3,2,1]);
        case 'YZ',
            h=permute(mean(values(:,idx,:),2),[3,1,2]);
        case 'XY',
            h=permute(mean(values(:,:,idx),3),[1,2,3]);
        case 'YX',
            h=permute(mean(values(:,:,idx),3),[2,1,3]);
    end

function h=plotMesh(o,varargin)
% h = grid.plotMesh(ax,[color],varargin) -- Plots wireframe of grid colored
% with var
%
% USAGE
% h = grid.plotMesh(ax,var,varargin);
% h = grid.plotMesh(var,varargin);
%
% h = handle to figure
% ax is handle to axes to plot onto.
% var is the variable which must be a 3D variable with data for all layers
%   and confining beds in sequence. If grid.LAYCBD is all zeros (no
%   confinging beds) then var could be HK VK or any layer variable.
%
% TO 120529

varargin(end+(1:2)) = {'visible','on'};

fsz = 12;
[ttl,varargin] = getProp(varargin,'title','mfilename unspecified variable');
[fsz,varargin] = getProp(varargin,'fonts' ,fsz);
[ax,varargin]  = getProp(varargin,'axis',[]);
if isempty(ax)
    [ax,varargin]  = getNext(varargin,'axis',gca);
end
[var,varargin] = getNext(varargin,'double',[]);

set(ax,'nextplot','add');

xlabel(ax,'x [m]','fontsize',fsz);
ylabel(ax,'y [m]','fontsize',fsz);
title(ax,ttl,'fontsize',fsz);

%% Now generate full size 3D arrays for all coordinates

if ~strmatchi('edgecolor',varargin), varargin(end+(1:2))={'edgecolor','k'};   end
if ~strmatchi('edgealpha',varargin), varargin(end+(1:2))={'edgealpha',0.25};   end

if ~isempty(var)
    
    if o.Nz ~= size(var,3)
        if size(var,3) == o.Nlay
            for i=o.Nlay:-1:1
                var(:,:,o.ITlay(i))=var(:,:,i);
            end
            for i=o.Ncbd:-1:1
                var(:,:,o.ITcbd(i))=NaN;
            end
        elseif size(var,3) == o.Ncbd
            for i=o.Ncbd:-1:1
                var(:,:,o.ITcbd(i))=var(:,:,i);
            end
            for i=o.Nlay:-1:1
                var(:,:,o.ITlay(i))=NaN;
            end
        else
            error('%s: size(var,3)=%d must equal gr.Nlay=%d or gr.Ncbd=%d',...
                mfilename,size(var,3),o.Nlay,o.Ncbd);
        end
    else
       % ok, use var as is
    end
    
    ix=1;      h(1) = surf(ax,squeeze(o.XGR(:,ix,:)),squeeze(o.YGR(:,ix,:)),squeeze(o.ZGR(:,ix,:)),squeeze(var(:,ix,:)));
    ix=o.Nx+1; h(2) = surf(ax,squeeze(o.XGR(:,ix,:)),squeeze(o.YGR(:,ix,:)),squeeze(o.ZGR(:,ix,:)),squeeze(var(:,max(1,ix-1),:)));
    iy=o.Ny+1; h(3) = surf(ax,squeeze(o.XGR(iy,:,:)),squeeze(o.YGR(iy,:,:)),squeeze(o.ZGR(iy,:,:)),squeeze(var(max(1,iy-1),:,:)));
    iy=1;      h(4) = surf(ax,squeeze(o.XGR(iy,:,:)),squeeze(o.YGR(iy,:,:)),squeeze(o.ZGR(iy,:,:)),squeeze(var(iy,:,:)));
    iz=o.Nlay+1; h(5) = surf(ax,squeeze(o.XGR(:,:,iz)),squeeze(o.YGR(:,:,iz)),squeeze(o.ZGR(:,:,iz)),squeeze(var(:,:,max(1,iz-1))));
    iz=1;      h(6) = surf(ax,squeeze(o.XGR(:,:,iz)),squeeze(o.YGR(:,:,iz)),squeeze(o.ZGR(:,:,iz)),squeeze(var(:,:,iz)));
else
    ix=1;      h(1) = surf(ax,squeeze(o.XGR(:,ix,:)),squeeze(o.YGR(:,ix,:)),squeeze(o.ZGR(:,ix,:)));
    ix=o.Nx+1; h(2) = surf(ax,squeeze(o.XGR(:,ix,:)),squeeze(o.YGR(:,ix,:)),squeeze(o.ZGR(:,ix,:)));
    iy=o.Ny+1; h(3) = surf(ax,squeeze(o.XGR(iy,:,:)),squeeze(o.YGR(iy,:,:)),squeeze(o.ZGR(iy,:,:)));
    iy=1;      h(4) = surf(ax,squeeze(o.XGR(iy,:,:)),squeeze(o.YGR(iy,:,:)),squeeze(o.ZGR(iy,:,:)));
    iz=o.Nlay+1; h(5) = surf(ax,squeeze(o.XGR(:,:,iz)),squeeze(o.YGR(:,:,iz)),squeeze(o.ZGR(:,:,iz)));
    iz=1;      h(6) = surf(ax,squeeze(o.XGR(:,:,iz)),squeeze(o.YGR(:,:,iz)),squeeze(o.ZGR(:,:,iz)));
end

if ~isempty(varargin)
    set(h,varargin{:});
end
drawnow;

%viewVector = get(ax,'CameraTarget') - get(ax,'CameraPosition');
%if viewVector(1)>0, setvisible(h([1 2])); else setvisible(h([2 1])); end
%if viewVector(2)>0, setvisible(h([3 4])); else setvisible(h([4 3])); end
%if viewVector(3)>0, setvisible(h([5 6])); else setvisible(h([6 5])); end

end

% function setvisible(h)
%     set(h(1),'visible','on');
%     set(h(2),'visible','off');
% end

function quiver(o,varargin)
%GRIDOBJ/QUIVER plots arrows on a grid
%
% USAGE:  gr.quiver([ax,]B,['power',power][,scale],plotoptions)
%
% power can be used to scale non linearly (it scale absolute values:
% F= sign(F) * |F|.^power.
%
% TO 131006

 [power,varargin]  = getProp(varargin,'power',1);
 [ax,varargin]     = getType(varargin,'axis',gca);
 [B ,varargin]     = getNext(varargin,'struct',[]);
 
 if isempty(B)
     error('%s: requires budget sctruct',mfilename);
 end
 

FR =  sum(B(end).term{strmatchi('flowR',B(end).label)},3); FR= FR./o.DY(:,:,1);
FF = -sum(B(end).term{strmatchi('flowF',B(end).label)},3); FF= FF./o.DX(:,:,1);
FR = [FR(:,1)  FR(:,1:end-1) FR(:,end-1)]; FR=(0.5*FR(:,1:end-1)+FR(:,2:end));
FF = [FF(1,:); FF(1:end-1,:);FF(end-1,:)]; FF=(0.5*FF(1:end-1,:)+FF(2:end,:));
if exist('power','var')
    FR = sign(FR).*abs(FR).^power;
    FF = sign(FF).*abs(FF).^power;
end

Idx = ~(FR==0 & FF==0);

quiver(ax,o.Xc(Idx),o.Yc(Idx),FR(Idx),FF(Idx),varargin{:});


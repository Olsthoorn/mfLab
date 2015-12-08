function plotBCN(gr,TYPE,BCN,iPer,iLay,lineSpec,varargin)
% gridObj.plotBCN(BCN,iPer,iLay,lineSpec,varargin);  % -% plot boundary condition
%   BCN is DRN RIV etc
%   iLay is layer (optional)
%   iPer is stress period number
%   simple minded approach, just for convenience.
%
% TO 121119

if nargin<5, lineSpec = 'r'; end
if nargin<4, iLay     = 1;   end
if nargin<3, iPer     = 1;   end
if nargin<2
    error('%s: must have at least one inputs [BCN/STRESS',mfilename);
end

if iscell(BCN)    
    BCN = cell2list(BCN(iPer));
else
    BCN = BCN(BCN(:,1)==iPer);
end

ixyz = BCN(BCN(:,2)==iLay,[4,3,2]);

if isempty(ixyz)
    msgId = 'mfLab:plotBCN:nothingToPlotInThisLayer';
    warning('on',msgId);
    warning(msgId,'%s: %s has no data for layer %d and stress period %d.',mfilename,TYPE,iLay,iPer);
    warning('off',msgId);
    return;
end

nRec= size(BCN,1);

x = NaN(6,nRec);
y = NaN(6,nRec);
z = NaN(6,nRec);
x(1,:) = gr.xGr(ixyz(:,1)'); x(2,:) = gr.xGr(ixyz(:,1)'+1);
x(3,:) = x(2,:); x(4,:) = x(1,:); x(5,:) = x(1,:);

y(1,:) = gr.yGr(ixyz(:,2)'); y(3,:) = gr.yGr(ixyz(:,2)'+1);
y(2,:) = y(1,:); y(4,:) = y(3,:); y(5,:) = y(1,:);

Idx = cellIndex(ixyz,gr.size);

z(1,:) = gr.ZM(Idx'); z(2,:)=z(1,:);
z(3,:) = z(1,:); z(4,:)=z(1,:); z(5,:)=z(1,:);

plot3(x(:),y(:),z(:),lineSpec,varargin{:});


function  [o] = changeLayer(o,layer,gr,IBOUND)
% This function will assign a new layer to the object. The layer needs to
% be defined. It will first assign the new layer and calculate
% corresponding values with it. Then the new values are interpolated to the
% their grid location. There is also a part which checks if the cells are
% within the ibound. This is turned off, because it happens automaticly.
%
% KG 171014

mesId = [class(o) ':pointsAtZeroIBOUND'];
warning('on',mesId);
            
for io = numel(o):-1:1
    jx   = o(io).iColx;
    jy   = o(io).iColy;
    [o(io).P, o(io).L] = gr.lineObjects( o(io).vertex(:,[jx jy]), layer);
    o(io).Idx = [o(io).P.idx];
                
    if nargin<2 || ~all( size(IBOUND(:,:,1))==o(io).grSize(1:2) ) || size(IBOUND,3)~=o(io).grSize(3)
        error('%s: need IBOUND in call and its size must be %s',...
            mfilename,sprintf(' %d',o(io).grSize));
    end

%     nEl1 = numel(o(io).Idx);
%     J = find(IBOUND(o(io).Idx)~=0);
%     o(io).cellLineVals    = o(io).cellLineVals(J,:);
%     o(io).P               = o(io).P(J);
%     o(io).A               = o(io).A(J);
%     o(io).sCell           = o(io).sCell(J);
%     
%     for i=1:numel(o(io).V)
%         if numel(o(io).V{i})>1
%             o(io).V{i}  = o(io).V{i}(J);
%         end
%     end
%     
%     o(io).Idx             = o(io).Idx(J);
%     nEl2 = numel(o(io).Idx);
%     
%     if nEl1~=nEl2
%         warning(mesId,'%s: %s type = %s, name = %s, numel cells %d > %d (%d vertices removed)',...
%         mfilename,class(o),o(io).type,o(io).name,nEl1,nEl2,nEl1-nEl2);
%     end
%     
%     if isempty(o(io).Idx)
%         warning(mesId,'%s: %s type = %s, name = %s is empty',...
%         mfilename,class(o),o(io).type,o(io).name);
%         warning(mesId,'off');
%         o(io) = [];
%     end
                
    %% Interpolate all values to their grid locations

    % cumulative line length at vertices
    line = o(io).vertex(:,[jx jy]);
    o(io).sLine = [0; cumsum(sqrt(sum(diff(line,1,1).^2,2)))];

    %Cumulative length along the intersected cells
    % to the mid of the intersection with the cell
    line = [[o(io).P.xm]; [o(io).P.ym]].';
    o(io).sCell = [0; cumsum(sqrt(sum(diff(line,1,1).^2,2)))];
    % elevation is known to P don't need to compute it

    % Interpolate the remaining values in the vertices to the cell centers
    % make sure the order of the values remains as they are in the
    % worksheet
    o(io).cellLineVals = interp1(o(io).sLine,o(io).vertex,o(io).sCell);

    %% Interpolate contour data to internal area
    [o(io).cellAreaVals,o(io).Idx] = ...
    gr.interpLaplacian([o(io).iColx,o(io).iColy],o(io).cellLineVals,[o(io).P.idx]);
    % add missing columns (id,x,y) using xm and ym
    o(io).A = gr.AREA3(o(io).Idx(:,1));
end     

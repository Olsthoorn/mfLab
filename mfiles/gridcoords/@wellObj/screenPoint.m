function o = screenPoint(o,gr,frac)
%% wellObj = wellObj.screenPoint(o,grSize,frac)
% Computes the z value, the layer number and the linear index number of a
% point on the screen at frac from its top.
% The results are in wellObj.UserData.screenPoint
%
% TO 120823

for i = 1:length(o)
    z = o(i).z(end)-frac*(o(i).z(end)-o(i).z(1)); %NA20130417

    Z = o(i).z(1)+[0; cumsum(o(i).DZ(:))]; %NA20130417

    ii = min(floor(interp1(Z,(1:length(Z))',z)),numel(o(i).DZ));
   

    o(i).UserData.screenPoint.z   = z;
    o(i).UserData.screenPoint.iLay= o(i).iLay(ii);
    o(i).UserData.screenPoint.LRC = o(i).LRC(ii,:);
    o(i).UserData.screenPoint.idx = o(i).idx(ii);
    o(i).UserData.screenPoint.idxMT3D = gr.Nxy*(o(i).LRC(ii,1)-1)+gr.Nx*(o(i).LRC(ii,2)-1)+o(i).LRC(ii,3);
end
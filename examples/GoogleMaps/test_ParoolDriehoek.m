% EXAMPLE showing how to get a google maps figure for given coordinates
%
% The coordinates of the LL and UR corner are specified in Dutch RD coords
%
% TO 110501 110514 121008

% Dutch national system coordinates
xGr = 122750 + (-1000:250:1000);
yGr = 485300 + (-1000:250:1000);
zGr = [0 -200];

gr = gridObj(xGr,yGr,zGr);

GM = googleMap(gr.xm([1 end]),gr.ym([1 end]),gr.zGr(end),'satellite','rd');

GM.image(xGr([1 end]),yGr([1 end]));


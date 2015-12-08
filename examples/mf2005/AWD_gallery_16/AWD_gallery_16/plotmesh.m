%plotmesh
%TO 110813
% script to plot the mesh
load underneath

figure; xlabel('x [m]'); ylabel('y [m]');
dphi=0.2; hrange=2:dphi:7; 
iLay=9;

grid.plot('xy','c-',1);
for i=1:length(waterbody)
    waterbody(i).plot('b-',2,'canal');
    waterbody(i).plot('r-',2,'pond');
    waterbody(i).plot('g-',1,'SpgPond');
end

grid.plot('xy','c-');
for i=1:length(waterbody),
    waterbody(i).plot('r',2,'canal');
    waterbody(i).plot('b',2,'pond');
    waterbody(i).plot('g',1,'SpgPond');
end
for i=1:length(Piez), Piez.plot('ro'); end
for i=1:length(Drain), Drain(i).plot('b',2); end

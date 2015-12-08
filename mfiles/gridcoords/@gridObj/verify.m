function verify(o,HK,tts)
%% gr.verify(tts,HK,Iy)  ---  Verify model by showing it
%% TO 120711 Trier

figure('position',(get(0,'screensize') + [50 50 0 0]) .* [1  1 0.75 0.75]);

xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');

view(3);

if ~exist('tts','var') || nargin<3,
    title('Verifying the Model');
else
    title(tts);
end

Iy = (1:o.Ny)';

h=o.plotMesh(log10(HK(Iy,:,:)),'facealpha',0.0); colormap;

set(h,'edgealpha',0.15);

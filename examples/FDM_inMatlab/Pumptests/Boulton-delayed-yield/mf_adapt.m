% is modflow perhaps more accurate?

AXIAL=1;

basename='fdmtest';

kh = 10;
sy = 0.1;
ss = 0.00001;
c  = 1;

%% Grid
xGr    = logspace(-2,5,71);
yGr    = [-0.5 0.5];
zGr    = ([-10:1:0 0.01:1.01])';

gr = gridObj(xGr,yGr,zGr,0,[],AXIAL);

HK     = gr.const(kh); kh(1,:,1)=0;
VKCB   = 1/c*ones(gr.Ny,gr.Nx,1); 
VK     = HK;
SY     = gr.const(sy);
SS     = gr.const(ss); SS(:,:,1)=sy./gr.DZ(:,:,1);

IBOUND = gr.const(1); IBOUND(:,:,1)=-1;
STRTHD = gr.const(0);

wel    = wellObj(basename,'wells',gr,HK,'PER');

T      = HK.*gr.DZ; T = T(1);
S      = sy;

save underneath T S
% is modflow perhaps more accurate?

AXIAL=1;

basename='fdmtest';

kh=10;
sy=0.1;
ss=0.1;

xGr=logspace(-2,5,71);
yGr=[-0.5 0.5];
zGr=[-1 0];

gr = gridObj(xGr,yGr,zGr);

HK = gr.const(kh);
VK = gr.const(kh);
SY = gr.const(sy);
SS = gr.const(ss);

IBOUND = gr.const(1); IBOUND(:,end,:)=-1;
STRTHD = gr.const(0);

well = wellObj(basename,'wells',gr,HK,'PER');

T   = HK.*gr.DZ;T=T(1);
S   = sy;

save underneath T S
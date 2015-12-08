% Seawat V4, Benchmark example, See Langevin e.a. 2008, p23ff

% Samani, N., M. Kompani-Zare, and D.A. Barry. 2004. MOD- FLOW
% equipped with a new method for the accurate simu- lation of
% axisymmetric flow. Advances in Water Resources 27, no. 1: 31?45.
%

%  TO 120118

basename='Samani';

% Assiging values from table for as far as needed for this model
kh     = 1e-5; % m/s
kv     = 1e-5; % m/s
Ss     = 1.03155e-3; % 1/m
strthd = 100.0; % m
%Q      = 6.28e-4; % m3/s
%R      = 11; % m
%R      = 39; % m

% if short
%     T = 19943; % seconds
%     tfac = 1.03663;
%     Nt=295;
%     t = factorspace(0.001,tfac,295);
% else
%     T=6.602e7;
%     tfac=1.58489;
%     Nt=36;
% end
%% Mesh using table data

AXIAL=1;

y1=8.0; f1=  NaN; % top of aquifer
y2=3.2; f2=1.387; % top of screen
y3=2.0; f3=1.400; % center of screen
y4=0.8; f4=1.410; % bottom of screen
y5=0.0; f5=  NaN; % bottom of aquifer

n=23; yT = factorspace(abs(y1-y2),f2,n);
n=17; yB =-factorspace(abs(y5-y4),f4,n);
n=18; yST=-factorspace(abs(y3-y2),f3,n);
n=18; ySB= factorspace(abs(y3-y4),f3,n);

zGr=unique([y1 y2 y3 y4 y5 y2+yT y4+yB y2+yST y4+ySB]); zGr=zGr(zGr>=0);

fprintf('%10.4f\n',min(diff(zGr)));

R=39; Nx=43; f=1.2350; xGr=factorspace(R,f,Nx); xGr(2)

gr = gridObj(xGr,[-0.5 0.5],zGr);

%% Generate all other matrices using table data

IBOUND   =gr.const(1);  IBOUND(:,end,:)=-1;
STRTHD   =gr.const(strthd);
HK       =gr.const(kh);
VK       =gr.const(kv);
SS       =gr.const(Ss);

%% Generate MNW objects and the acompanying PNTSRC for MT3DMS

MNW = MNW1Obj(basename,'wells',gr,HK);


save Underneath Ss
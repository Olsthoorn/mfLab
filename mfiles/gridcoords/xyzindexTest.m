function xyzindexTest(usage)
%XYZINDEXTEST tests function xyzindex
%    xyzindex computes the cell indices of points in 1D, 2D or 3D grid
%
% Example:
%    xyzindexTest(usage); % usage 1..8
%
% Usage is an integer between 1 and 8, used in switch to test different cases
%
% todo: test with LAYCBD
%
% SEE ALSO: also IdxMatlab2Modflow cellIndex xyzindex linegrid inpolygon
%
% TO130127

xyzIsRandom = false;

if nargin<1 || usage<1 || usage>9
    error('%s: usage must be between 1 and 8, not %d',usage);
end

zT  =    0;
zB  = -100;

xGr = 0:100:1000;          Nx = length(xGr)-1;
yGr = (1000:-100:0)';      Ny = length(yGr)-1;

zTop = zT * ones(Ny,Nx); %zTop = zTop + 25*(rand(size(zTop))-0.5);
zBot = zB * ones(Ny,Nx); %zBot = zBot + 25*(rand(size(zBot))-0.5);

Nz = 10;
Z = NaN(Ny,Nx,Nz+1); Z(:,:,1)=zTop; Z(:,:,end)=zBot;
for iy=1:Ny
    for ix=1:Nx
        Z(iy,ix,:) = interp1([0 Nz],[zTop(iy,ix) zBot(iy,ix)],0:Nz);
    end
end

gr = gridObj(xGr,yGr,Z);

if xyzIsRandom
    xyz = rand(20,3);
else
    xyz = [
        0.7092    0.0990    0.9554
        0.6413    0.5710    0.1427
        0.1741    0.3259    0.5126
        0.0622    0.4505    0.9719
        0.4067    0.5778    0.6483
        0.4631    0.0748    0.6147
        0.2027    0.0573    0.4697
        0.8695    0.3010    0.5778
        0.5979    0.5217    0.9113
        0.0230    0.5619    0.3762
        0.8994    0.2416    0.2288
        0.4529    0.9127    0.4235
        0.0580    0.8257    0.2736
        0.1063    0.4445    0.4446
        0.9984    0.9821    0.6275
        0.8663    0.5783    0.5346
        0.6152    0.2344    0.3854
        0.0269    0.8106    0.8735
        0.3225    0.4513    0.3003
        0.4638    0.2500    0.4000];
end

xyz = xyz .* (ones(20,1)*[xGr(end),yGr(1),Z(end)]);

x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);

% USAGES:
switch usage
    case 1
       [ix,xR]             = xyzindex(x,xGr);
       display([x ix xR]);
    case 2
       [ix,iy,xR,yR]       = xyzindex(x,y,xGr,yGr);
       display([x ix xR y iy yR]);
    case 3
       [ix,iy,xR,yR]       = xyzindex([x y],xGr,yGr);
       display([x ix xR y iy yR]);
    case 4
       [ix,iy,iz,xR,yR,zR] = xyzindex(x,y,z,xGr,yGr,Z);
       display([x ix xR y iy yR z iz zR]);
    case 5
       [ix,iy,iz,xR,yR,zR] = xyzindex([x,y,z],gr);
       display([x ix xR y iy yR z iz zR]);
    case 6
         Idx = xyzindex([x,y,z],gr);
         display(Idx);
    case 7
         Idx = xyzindex([x,y,z],xGr,yGr,Z);
         display(Idx);
    case 8
         Idx = xyzindex(x,y,z,xGr,yGr,Z);
         display(Idx);
    case 9
        LAYCBD = [1 0 0 1 1 1 0 1];
        gr = gridObj(xGr,yGr,Z,LAYCBD);
       [ix,iy,iz,xR,yR,zR] = xyzindex([x,y,z],gr);
       display([x ix xR y iy yR z iz zR]);
end

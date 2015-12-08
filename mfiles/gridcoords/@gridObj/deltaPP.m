function dPP = deltaPP(o,well,KR,KV)
%DELTAPP analytically computes extra drawdown due to partial penetration given well and anquifer
%
% Example:
%     dPP = grid.deltaPP(well,HK[,VK])
% 
%     The grid o must be axially symmetric (o.AXIAL<>0);
%     The grid must have uniform layers.
%     If o.Ny>1, there must be a well in the first column of each row only.
%     Then make sure CHANI=1e-8 to prevent shortcuts between rows (your
%     responsibility)
%
%     if VK is specified K = (HK.*HK.*VK)^(1/3) is used.
%
% See also: hantush hantushn hantushE
%
% TO 120930

if ~o.AXIAL, error('%s: Grid must be axially symmetric',mfilename); end
if ~o.layersAreUniform,
    error('%s: layers in grid must be uniform');
end

if any([well.ix]~=1), error('%s: all wells must be in col 1',mfilename); end
if numel(well) ~= o.Ny
    error('%s: the number of wells must match the number of rows in the grid',mfilename);
end
if ~all( (1:o.Ny) == [well.iy] )
    error('%s, there must be one well per row',mfilename);
end
    
%% Verify steady drawdown by subtracting the drawdown at r<5000 m for t = 100 days wit the analytical drawdown
kr = mean(mean(KR,3),2);

if nargin>3,
    kv = mean(mean(KV,3),2);
    k = (kr.*kr.*kv).^(1/3);
    rStretchFac = sqrt(k./kr);
    zStretchFac = sqrt(k./kv);
else
    k = kr;
    rStretchFac = ones(size(k));
    zStretchFac = ones(size(k));
end

D      = sum(o.dz) *zStretchFac;

for iy     = o.Ny:-1:1
    Q(iy)  = well(iy).Q;
    a(iy)  = abs(well(iy).z(  1)) * zStretchFac(iy);
    b(iy)  = abs(well(iy).z(end)) * zStretchFac(iy);
    d(iy)  = abs(a(iy)-b(iy));
    dq(iy) = Q(iy)/(pi*k(iy)) /(pi*d(iy));
end

XM = o.XM; for iy=1:o.Ny, XM(iy,:,:)=XM(iy,:,:)*rStretchFac(iy); end
ZM = o.ZM; for iy=1:o.Ny, ZM(iy,:,:)=ZM(iy,:,:)*zStretchFac(iy); end

dPPOld = o.const(0);
dPP    = o.const(0);
for iy=1:o.Ny
    for n=1:250
        dPP(iy,:,:) = dPP(iy,:,:) + dq(iy) * ...
            1/n*(sin(n*pi*b(iy)/D(iy)) - sin(n*pi*a(iy)/D(iy))).* ...
            cos(n*pi*ZM(iy,:,:)/D(iy)).*besselk(0, n*pi*XM(iy,:,:)/D(iy));
        if rem(n,50)==0
            fprintf('iy=%2d, iteration %2d: %g\n',iy,n,max(max(max(abs(dPP-dPPOld)))));
        end
        dPPOld(iy,:,:) = dPP(iy,:,:);
    end
    fprintf('\n');
end


function [dtadv,dtdisp,dtss,Ladv,Ldisp,Ltss] = DtRequirements(gr,basename,PEFF,it)
% [Dtadv,Dtdisp,Dtss] = DtRequirements(basename,PEFF);
%
% TO 120904

if nargin<4, it=1; end

if gr.Nx<2 || gr.Ny<2 || gr.Nz<2
    warning('%s: only works if all dimensions are >=2'); %#ok
    dtadv   = NaN;
    dtdisp  = NaN;
    dtss    = NaN;
    Ladv    = NaN;
    Ldisp   = NaN;
    Ltss    = NaN;
    return;
end

% get budget file with cell2cell flows
B = readBud([basename '.BGT']);

% get MT3D parameters
[mt3dHdr,mt3dvals] = getExcelData(basename,'MT3D','vertical');

% get Layer parameters
[layHdr ,layVals ] = getExcelData(basename,'LAY'  ,'horizontal');

% Get the largest diffusion coefficient in the model over all compoments
NCOMP = mt3dvals(strmatchi('NCOMP',mt3dHdr));
ID    = strmatchi('DMCOEF',layHdr);
D     = max(max(layVals(:,ID(1:NCOMP))));

% flow across cel faces
FR = B(it).term{strmatchi('FLOWRIGHTFACE',B(it).label)};
FL = B(it).term{strmatchi('FLOWLOWERFACE',B(it).label)};
FF = B(it).term{strmatchi('FLOWFRONTFACE',B(it).label)};

% absolute maximum velocity in all cells
vx = 0.5 * ( abs(cat(2,FR(:,1,:),FR(:,1:end-1,:))) + abs(cat(2,FR(:,1:end-1,:),FR(:,end-1,:))) );
vy = 0.5 * ( abs(cat(1,FF(1,:,:),FF(1:end-1,:,:))) + abs(cat(1,FF(1:end-1,:,:),FF(end-1,:,:))) );
vz = 0.5 * ( abs(cat(3,FL(:,:,1),FL(:,:,1:end-1))) + abs(cat(3,FL(:,:,1:end-1),FL(:,:,end-1))) );

% velocities
vx = vx./(gr.DY.*gr.DZ.*PEFF);
vy = vy./(gr.DX.*gr.DZ.*PEFF);
vz = vz./(gr.DX.*gr.DY.*PEFF);

% assume max velocity in any cell
vabs = sqrt(vx.^2+vy.^2+vz.^2);

% dispersion coefficients
DL = D + vabs.*gr.const(layVals(:,strmatchi('AL',layHdr,'exact')));
DT = D./strmatchi('TRPT',layHdr,'exact');
DV = D./strmatchi('TRPV',layHdr,'exact');

% requirements in every cell
Dtadv  = 1./(vx./gr.DX + vy./gr.DY + vz./gr.DZ);
Dtdisp = 0.5 ./ ( (DL+DT+DV).*(1/gr.DX.^2 + 1./gr.DY.^2 + 1./gr.DZ.^2));

I = find(~ismember(B(it).label,{'CONSTANTHEAD','FLOWRIGHTFACE','FLOWLOWERFACE','FLOWFRONTFACE'}));

Qcell = abs(B(it).term{strmatchi('CONSTANTHEAD',B(it).label)});
if ~isempty(I)
    for i=1:length(I)
        Qcell = Qcell + abs(B(it).term{I(i)});
    end
end

Dtss   = gr.Vlay .* PEFF./Qcell;

% overall requirements
dtadv  = min(Dtadv(:));        Idx1 = find(Dtadv == dtadv);  Idx1=Idx1(1);
dtdisp = min(Dtdisp(:));       Idx2 = find(Dtdisp== dtdisp); Idx2=Idx2(1);
dtss   = min(Dtss(Dtss(:)>0)); Idx3 = find(Dtss  == dtss);   Idx3=Idx3(1);

Ladv  = cellIndices(Idx1,gr.size,'LRC');
Ldisp = cellIndices(Idx2,gr.size,'LRC');
Ltss  = cellIndices(Idx3,gr.size,'LRC');


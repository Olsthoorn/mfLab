function PNTSRC = getPNTSRC(BCN,type,STCONC)
%% PNTSRC = getPNTSRC(BCN,type,STCONC) -- generate PNTSRC for this boundary condition using STCONC
% BCN is boundary condition array like WEL, DRN, RIV GHB, CHD having
% [SP L R C ] as the first 4 columns
% type is a 3 letter string like 'WEL', 'DRN', GHB' etc
% STCONC is either a full 3D numerical array or a cell array where every
% cell contains a full 3D concentration array.
% any array of concentration may be used al long as the size matches
% [Ny Nx Nz]
%
% TO 120511

if isempty(BCN)
    PNTSRC=[];
    return;
end

Idx = cellIndex(BCN(:,[4 3 2]),size(STCONC));
Idx =Idx(~isnan((Idx)));

switch upper(type)
    case 'CCC', ITYPE=-1; % constant concentration cell
    case 'CHD', ITYPE= 1; % constant head cell
    case 'WEL', ITYPE= 2; % well
    case 'DRN', ITYPE= 3; % drain
    case 'RIV', ITYPE= 4; % river
    case 'GHB', ITYPE= 5; % general head boundary
    case 'RCH', ITYPE= 7; % recharge
    case' EVT', ITYPE= 8; % evapotranspiration
    case 'MLS', ITYPE=15; % mass loading source       
    case 'STR', ITYPE=21; % stream routing
    case 'RES', ITYPE=22; % reservoir
    case 'SFH', ITYPE=23; % Specified flow and head boundary
    case 'IBS', ITYPE=24; % Inter-bed storage
    case 'TRP', ITYPE=25; % Transient leakage
    case 'LAK', ITYPE=26; % Lake
    case 'MNW', ITYPE=27; % Multi-node well
    case 'DRT', ITYPE=28; % Drain with return flow
    case 'ETS', ITYPE=29; % Segmented evapotranspiration
    case 'HSS', ITYPE=50; % HSS Mass Loading
    otherwise
        error('mfLab:getPNTSRC:unknownType',...
            'Unknown or reserved type=%s',type);
end

if isempty(Idx)
    PNTSRC={[]};
else
    if iscell(STCONC)
        C=NaN(size(Idx,1),length(STCONC));
        for j=1:size(STCONC)
            C(:,j)=STCONC{j}(Idx);
        end
        PNTSRC = {[BCN(:,1:4) C(:,1) ones(size(Idx))*ITYPE C]};
    else
        PNTSRC = {[BCN(:,1:4) STCONC(Idx) ones(size(Idx))*ITYPE]};
    end
end
function List=BCN_Struct(s,S)
%BCN_STRUCT makes a list of form [iPer iLay iRow iCol rest]
%
% Example:
%    List=BCN_Struct(s,S)
%
% INPUT:
%   rest depends on the string S, which may be one of 'WEL' 'DRN' 'RIV' 'GHB' ...
%
%   The function is used to generate boundary condition (stress) lists from the
%   struct s. It is used to construct the groundwater model of
%   the Amsterdam Water Supply
%
% Used in:
%   examples/swt_v4/AWD/AMWADU/amwadu2mfLab/AMWADUdata
%
% TO 100528

tic
fprintf('Busy extracting %s ... ',S);

NPER=length(s.itmp);
NCEL=length(s.row);

switch S
    case 'WEL', L=7; % [iPer iL iR iC Q Qtot NofWells] 
    case 'DRN', L=6; % [iPer iL iR iC Elev Cond]
    case 'DRT', L=6; % [iPer iL iR iC Elev Cond]
    case 'GHB', L=6; % [iPer iL iR iC Head Cond] 
    case 'RIV', L=7; % [iPer iL iR iC Head Cond RBot]
    otherwise, error('Illegal arguments %s in BCN_struct',S);
end

%%Prevent dynamic memory allocation, first count lenght of list

List=NaN(NPER*NCEL,L);

for iPer=1:NPER
    J=(iPer-1)*NCEL;
    List(J+(1:NCEL),1)  = iPer;
    List(J+(1:NCEL),3)  = s.row;
    List(J+(1:NCEL),4)  = s.col;
    
    switch S
        case 'WEL', K=find(~isnan(s.flow(:,iPer)))';
        otherwise,  K=find(~isnan(s.head(:,iPer)))';
    end

    for iO=K
        I=J+s.idnr{iO};
        
        List(I,2)=s.layer(iO);
        switch S
            case 'WEL'
                List(I,5)=s.flow(iO,iPer)./s.wells(iO);
                List(I,6)=s.wells(iO);
                List(I,7)=s.flow(iO,iPer); % Total flow for total object
                List(I,8)=s.gis_id(iO);
            case 'DRN'
                if isfield(s,'length'),
                    List(I,5)=s.head(iO,iPer);
                    List(I,6)=s.length(iO)./s.resistance(iO);
                elseif isfield(s,'mp')
                    List(I,5)=s.mp(iO);
                    List(I,6)=s.area(iO)./s.resistance(iO);
                else
                    List(I,5)=s.head(iO,iPer);
                    List(I,6)=s.area(iO)./s.resistance(iO);
                end
                List(I,7)=s.gis_id(iO);
            case 'GHB'
                List(I,5)=s.head(iO,iPer);
                List(I,6)=s.area(iO)/s.resistance(iO);
                List(I,7)=s.gis_id(iO);
            case 'RIV'
                List(I,5)=s.head(iO,iPer);
                List(I,6)=s.area(iO)/s.resistance(iO);
                List(I,7)=s.depth(iO);
                List(I,8)=s.gis_id(iO);
        end
    end
end

List=List(~isnan(List(:,5)),:);

fprintf('%10d %s in %10g seconds\n',size(List,1),S,toc);

function OBJ=BCN_Curb(OBJ,ix,iy,iz)
%BCN_CURB cuts out and renumber BCN when cutting out peice of larger model
%
% Example:
%    OBJ=BCN_Curb(OBJ,ix,iy,iz)
%
% Make sure that a boundary condition list of the form
% [iPer iLay iRow iCol ....]
%
% Used to generate a child model from the model of the Amsterdam Water
% Supply
%
% TO 100528

ix=[min(ix),   max(ix)];
iy=[min(iy),   max(iy)];
iz=[min(iz(:)),max(iz(:))];

OBJ(OBJ(:,2)<iz(1) | OBJ(:,2)>iz(2),1)=NaN;
OBJ(OBJ(:,3)<iy(1) | OBJ(:,3)>iy(2),1)=NaN;
OBJ(OBJ(:,4)<ix(1) | OBJ(:,4)>ix(2),1)=NaN;

OBJ=OBJ(~isnan(OBJ(:,1)),:);

OBJ(:,2)=OBJ(:,2)-iz(1)+1;
OBJ(:,3)=OBJ(:,3)-iy(1)+1;
OBJ(:,4)=OBJ(:,4)-ix(1)+1;

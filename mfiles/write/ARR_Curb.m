function OBJ=ARR_Curb(OBJ,iy,ix,iz)
% OBJ=ARR_Curb(OBJ,iy,ix,iz)
% Make sure that a model Array of the form [Ny,Nx,Nz]
% is reduced to the limites of the child model given by
% ix, iy and iz indices in the parent model
% Used to generate a child model from the model of the Amsterdam Water
% Supply
%
% TO 100528
   ix=[min(ix),   max(ix)];
   iy=[min(iy),   max(iy)];
   iz=[min(iz(:)),max(iz(:))];
   
   OBJ=OBJ(iy(1):iy(2),ix(1):ix(2),iz(1):iz(2));

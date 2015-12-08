function [isLay,LAYCBD]=isLayer(Nz,LAYCBD)
%ISLAYER generates pointers to actual model node layers in the Z-implied layers
% It also completes the vector LAYCBD needed by MODFLOW to discern model
% node layers from confining beds.
%
%  USAGE:  [isLay, LAYCBD] = isLayer(Nz,LAYCBD)
%
%  Nz is the number of layers implied in the Z-array with which gridObj was
%  called.
%  LAYCBD is a vector Nlay long with a value > 0 for every layer that
%         has a confing bed below it. It may be incomplete, e.g. a single
%         zero. It will be completed in this function as needed.
%
%  Used in gridObj
%
%  See also: gridObj
%
% TO 120414 151115

% How many model node layers and confining beds do we have?
LAYCBD = logical(LAYCBD(:)) .'; % hor vector
Ncbd   = sum(LAYCBD);
Nlay    = Nz - Ncbd;

%% Clean up LAYCBD
LAYCBD(Nlay+1:end)=[];          % truncate when LAYCBD is too long
LAYCBD(end+1:Nlay)=false;       % complete when LAYCBD is too short
LAYCBD(end)       =false;       % last LAYCBD must always be zero

% Manipulate to get the pointers right
zImplied = [ones(size(LAYCBD)); LAYCBD]; % complete with model layers as first column
zImplied = cumsum( zImplied(:));            % count z-implied layers
zImplied = reshape(zImplied,[2,Nlay]);
zImplied = zImplied(1,:)';

if Nz ~= max(zImplied(:))
    error('Something wrong here, Nz=%d not equal to %d !',Nz,max(zImplied(:)));
end

isLay                = false(Nz,1);
isLay(zImplied) = true; % true of z-layer is a model layer

LAYCBD = LAYCBD .';



function [L2N, C2N] = convert2Natural(o)
%% [L2N, C2N] = convert2Natural(LAYCBDold,LAYCBDnew)
%
% Convert iLay or iCbd to iz the natural stack order
%
% L2N  -> layer sequence to natural sequence
% C2N  -> cbd   sequence to natural seqence
%
% iz = L2N(iLay);

% TO 120507

% Split LAYCBD vector into 2 columns, col1 representing the layers and col1 the cbd or 0 if absent

%% From input side generate the 2-column array with natural stack numbers of the input layers and cbd

% in two steps: first split and then put the numbers, leaving zeros where cbd absent
Ifrom = [ones(size(o.LAYCBD')); o.LAYCBD'];  % i.e.: [1 0 1 1 0 0] --> [1 1 1 1 1 1; 1 0 1 1 0 0]
Ifrom = reshape( (Ifrom(:)>0) .* cumsum(Ifrom(:)),size(Ifrom))'; %         --> [1 3 4 6 8 9; 2 0 5 7 0 0]'

% natural stack numbers of model layers
L2N = Ifrom(:,1);

% natural stack numbers of model cbd
C2N= Ifrom(Ifrom(:,2)>0,2);

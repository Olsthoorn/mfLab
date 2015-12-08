function [L2N, C2N, LC2N, N2L, N2C] = convertLAYCBD(LAYCBDold,LAYCBDnew)
%CONVERTLAYCBD converts layers from LAYCBDold to LAYCBDnew (with selftest)
%
% Example:
% [L2N, C2N, LC2N, N2L, N3C] = convertLAYCBD(LAYCBDold,LAYCBDnew)
%
% Generates indices of model layers and confining beds between two grids
% onw with LAYCBSold and one with LAYCBDnew. It is used in methods of
% modelObj to convert one MODFLOW grid into another having different
% distributions of model layers and confining beds. See further narrative
% below.
%
% L2N  -> layer sequence to natural sequence
% C2N  -> cbd   sequence to natural seqence
%
% LC2N -> for values like HK having no CBD, the input sequence to natural sequence
%  e.g. the natural layer corresponging to cbd 2 uses HK of layer above this CBD
%
% N2L  -> natural sequence to output layer sequence
% N2C  -> natural sequence to output cbd   sequence
%
% convertLAYCBD for selftest
%
% Get index to convert from old layer to new layers given
% LAYCBDold and LAYCBDnew
%
% Also compute the location in the new layer stack in case the old layer
% array and the old cbd array are given separately.
%
% See also: modelObj
%
% TO 120507

selftest = nargin==0;

if selftest
  LAYin = [11 13 14 16 18 19 21 22 23 25 27]' % selftest input values for layers
  CBDin = [12 15 17 20 24 26]'                % selftest input values for CBD
  LAYCBDold = [1 0 1 1 0 1 0 0 1 1 0]'        % test input LAYCBD
  LAYCBDnew = [0 0 1 0 1 1 0 0 1 1 0 0]'      % test output LAYCBD
end

% How many old and new layers and cbd do we have (use LAYCBDs to find out)
NlayOld = length(LAYCBDold);  NlayNew = length(LAYCBDnew);
NcbdOld = sum(LAYCBDold>0);   NcbdNew = sum(LAYCBDnew>0);

%% Check if sum of layers and cbd remains constant
if NlayOld+NcbdOld ~= NlayNew+NcbdNew
    error(['NlayOld(%d)+NcbdOld(%d) ~= NlayNew(%d)+NcbdNew(%d)\n',...
           'where Nlay=length(LAYCBD) and Ncb=sum(LAYCBD>0)'],...
           NlayOld,NcbdOld,NlayNew,NcbdNew);
end

%% Natural stack or layer order
% The natural stack or layer order is the one for the case without any
% confining beds, that is, when all all(LAYCBD==0).
%
% To transfer old layers and old cbd to new layers and new cbd we
% attributed them natural stack numbers, giving a unique layer number to
% each of the old and the new layers.
%
% Having done this, we can immediately place old layers and cbd in their
% natural stack order. And also we can immediately thereafter pick the new
% layers and cbd from the natural stack.
%
% However, with each intput and output layer and cbd having its unique
% natural stack number, we can also immediately transfer input layers and
% cbd to their correct output layers and cbd, circumventing the generation
% of a natural layer stack altogether.
%
% Transfer from the input  to the natural layer stack
% Transfer from the output to the natural layer stack
% When done put natural stack order back to back in sequence to transfer
% from input layer sequence to output layer sequence. Note that in case we
% have CBD on eithe side in arbitrary numbers and order we need the natural
% stack as intermediate. This is what we, therefore, always apply.
%

%% Input side to natural layer stack
% Ioriginal is the natural stack numbers attributed to the layers and cbd
% of the input, which we combine as [LAYold ; CBDold]. Ioriginial is
% implied by LAYCBDold.

%% First step
% All change LAYCBD into and cumsum, while afterwards set values zero where LAYCBD zero.
% Try this

% The first idea is to Split LAYCBD vector into 2 columns, first representing the layers and second the cbd or 0 of absent

%% From input side generate the 2-column array with natural stack numbers of the input layers and cbd
% in two steps: first split and then put the numbers, leaving zeros where cbd absent
Ifrom = [ones(size(LAYCBDold(:)')); LAYCBDold(:)'];  % i.e.: [1 0 1 1 0 0] --> [1 1 1 1 1 1; 1 0 1 1 0 0]
Ifrom = reshape( (Ifrom(:)>0) .* cumsum(Ifrom(:)),size(Ifrom))'; %         --> [1 3 4 6 8 9; 2 0 5 7 0 0]'

% 1st col input sequence (1:Nlay), 2nd col natural stack numbers of the input layers
L2N = [(1:NlayOld)', Ifrom(:,1)];         % Source layer -> natural     --> [1 2 3 4 5 6; 1 3 4 6 8 9]'

% 1st col input sequence (1:Ncbd), 2nd col natural stack numbers of the input cbd
C2N= [(1:NcbdOld)' Ifrom(Ifrom(:,2)>0,2)];  % Sourc  CBD -> natural     --> [1 2 3 ; 2 5 7];

% 1st col input sequence of layer list, 2nd column natural sequence of CBND
LC2N =Ifrom(Ifrom(:,2)>0,:); % to use LAY values for CBD that don't have them (PEFF, STRTHD whatever)

%% From output side the same procedure
Ito   = [ones(size(LAYCBDnew(:)')); LAYCBDnew(:)'];
Ito   = reshape( (Ito(:)>0) .* cumsum(Ito(:)),[2,NlayNew])';

% 1st col natural stack nrs, 2nd col output sequence (1:Nlay)
N2L   = [Ito(:         ,1) (1:NlayNew)'];    % Target layer --> natural

% 1st cal natural stack nrs, 2nd col output sequence (1:Ncbd)
N2C   = [Ito(Ito(:,2)>0,2) (1:NcbdNew)'];    % Target CBD   --> natural

%% Do selftest if nargin==0
if selftest

    Natural = NaN(NcbdOld+NlayOld,1)      %#ok<*NOPRT> % natural stack allocation

    LAYout = NaN(NlayNew,1);              % output allocation
    CBDout = NaN(NcbdNew,1);              % output allocation

    Natural(L2N( :,2)) = LAYin(L2N( :,1)) % input layer values to natural stack
    Natural(C2N(:,2)) = CBDin(C2N(:,1)) % input cbd   values to natural stack
  % Natural(LC2N(:,2)) = LAYin(LC2N(:,1)) % remplacent CBD values if no actual values exist in CBP
 
    LAYout(N2L(:,2)) = Natural(N2L(:,1)) %#ok<NASGU> % output layer values from natural stack
    CBDout(N2C(:,2)) = Natural(N2C(:,1)) %#ok<NASGU> % output cbd   values from natural stack
end
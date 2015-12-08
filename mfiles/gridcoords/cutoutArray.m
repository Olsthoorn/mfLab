function [Alay,Acbd] = cutoutArray(A_In,Ix,Iy,LAYCBDold,LAYCBDnew)
%CUTOUTARRAY cutout 3D subarray using indices Ix,Iy and definition of confing beds of input and output
%
% Example:
%    [AlayOut,AcbdOut] = cutoutArray(A_In,Ix,Iy,LAYCBDold,LAYCBDnew)
%
% OUTPUTS:
%     Alay   Portion of 3D array (only model layers) based on the
%            selection Ix,Iy and LAYCBDold and LAYCBDnew
%     Acbd   Same, but only the confining beds.
%
% Used in conjuntion with modelObj to refine and coarsen model grids
%
% INPUTS:
% A_in:        is the aquifer array with aquitard array concatenated into the 3rd
% dimension:  for instance A_In = cat(HK,VKCB,3); or just the aquifer data
%             if in the new model the aquitard data can be used for aquifer data.
%             Then A_in = HK suffices, like with A_In = PEFF or A_In= Ss.
%
%    If AcbdIn is not specified the values of AlayOut will be used.
%    Old layer nrs are deduced from AlayOut.
%    The total number of layers and aquitards is deduced from
%    length(LAYCBD)+sum(LAYCBD>0), where LAYCBD reflects the new situation
%    to reflect the new coordinate indices Ix and Iy
%
% These function are superseded by methods of the gridObj
%
% see also: gridObj modelObj cutoutBCN cutoutXTA RefineBCN RefindGrid removeCBD
%
% TO 120507

%% Check input
if min(Ix)<1|| max(Ix)>size(A_In,2),
    error('%s Ix indices must be >0 and <%d',mfilename,size(A_In,2)); end

if min(Iy)<1 || max(Iy)>size(A_In,1),
    error('%s Iy indices must be <0 and <%d',mfilename,size(A_In,1)); end

NlayOld = length(LAYCBDold);
NcbdOld = sum(LAYCBDold>0);

NlayNew = length(LAYCBDnew);
NcbdNew = sum(LAYCBDnew>0);

if NlayOld+NcbdOld ~= NlayNew+NcbdNew
    error(['NlayOld(%d)+NcbdOld(%d) ~= NlayNew(%d)+NcbdNew(%d)\n',...
           'where Nlay=length(LAYCBD) and Ncb=sum(LAYCBD>0)'],...
           NlayOld,NcbdOld,NlayNew,NcbdNew);
end

if  ~(size(A_In,3) == NlayOld+NcbdOld || size(A_In,3) == NlayOld),
    error(['%s sum of layers in the input array must be equal length(LAYCBD)=%d or\n',...
        'equal to length(LAYCBD)+sum(LAYCBD>0)=%d.'],...
        mfilename,NlayOld,NlayOld+NlayNew);
end

%% Input Side

Ioriginal= layerOrder(LAYCBDold,NlayOld);

if size(A_In,3)< NlayOld + NcbdOld
    Ioriginal(NlayOld+(1:NcbdOld))=Ioriginal(NlayOld+(1:NcbdOld))-1;
end
    
I    = sortrows([Ioriginal  (1:NlayOld+NcbdOld)']);
    
if size(A_In,3)< NlayOld+ NcbdOld
    for i=1:size(I,1)-1; if I(i+1,1)==I(i,1), I(i+1,2)=I(i,2); end; end
end

I = I(:,2);   % Iin(i) is the layer in the orginal stack of layer(i) in input

%% Output side


Ioriginal= layerOrder(LAYCBDnew,NlayNew);
    
J    = sortrows([Ioriginal,   (1:NlayNew+NcbdNew)']);
    
J = J(:,2);   % Iin(i) is the layer in the orginal stack of layer(i) in input

Jout = sortrows([J I]);

Jout = Jout(:,2);

%% finalize

Alay = A_In(Iy,Ix,Jout(1:NlayNew));
if nargout>1,
    Acbd = A_In(Iy,Ix,Jout(end-NcbdNew+(1:NcbdNew)));
end

end

function Io = layerOrder(LAYCBD,Nlay) 
    Io = [ones(size(LAYCBD(:)')); LAYCBD(:)'];
    Io = reshape(Io(:) .* cumsum(Io(:)),[2,Nlay])';
    Io = [Io(:,1); Io(Io(:,2)~=0, 2)];
end


function [isaquif,LAYCBD]=isAquifer(NZ,LAYCBD)
%ISAQUIFER logical vector telling which model layer is an aquifer and which is not
%
% Example:
%    [isaquif,LAYCBD] = isAquifer(NZ,LAYCBD)
%
% where
%    isaquif (Nlay+Ncbd) is a vector with 1 if layer iand 0 if confining bed
%
% used by: gridObj
%
% Notice:
%    LAYCBD (Nz,1) is 0 if layer has no confining bed and 1
%       is has one beneath it.
%    NZ = totale number of layers, aquitard + aquifers
%
% See also: gridObj
%
% TO 120406

laycbd=LAYCBD(:); 
Nz1=length(laycbd)+sum(laycbd>0);

if Nz1<NZ, laycbd=[laycbd; zeros(Nz1-NZ,1)]; end

laycbd(end)=0;
    
while 1
    Nz1   =  length(laycbd)+sum(laycbd>0);
    if Nz1>NZ,
        laycbd(end)=[];
        laycbd(end)=0;
    else
        LAYCBD=laycbd;
        break;
    end
end

% LAYCBD is already compatible here, so no checking necessary
isaquif=ones(length(LAYCBD)+sum(LAYCBD>0),1);      % layer indices
k=0;             % counter
for iz=1:length(LAYCBD) % is number of aquifer layers or node containing layers
    k=k+1;
    if LAYCBD(iz)>0
        k=k+1;
        isaquif(k)=0;
    end
end

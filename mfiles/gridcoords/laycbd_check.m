function LAYCBD=laycbd_check(Nz,LAYCBD)
%LAYCBD_CHECK checks consistency if LAYCBD with number of layers in MODFLOW model
%
% Example:
%    LAYCBD=laycbd_check(Nz,LAYCBD)  
%
% Nz is number of layers including confining beds (size(zm,3) or size(zGr,3)-1
%
% See also: gridObj
%
% TO 120407

if Nz<length(LAYCBD)+sum(LAYCBD>0), LAYCBD(end)=0; end

if Nz~=length(LAYCBD)+sum(LAYCBD>0)
    
    error('mfLab:laycbd_check:sizeLAYCBD',...
        ['Nz= %d (aquifers+confining beds) does not match LAYCBD\n',...
         'as total number of aquifers=%d + aquitards=%d must equal,Nz=%d\n'],...
         Nz,length(LAYCBD),sum(LAYCBD>0),Nz);
end

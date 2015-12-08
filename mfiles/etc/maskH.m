function H=maskH(H)
%MASKH masks heads H(i).values using HINACT.mat saved by mf_setup
%
% USAGE:
%    H = maskH(H)
%
%    where H is the struct provided by readDat([basename '.HDS'];
%
% See Also  maskBud readDat readMT3D readBud maskHC maskC
%
% TO 090315 100420 100718 100911 120409 120509 120524

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin<1, error('mflab:maskH:insufficientInputArguments',...
    'maskH: 1 input argument required\n');
end

try
    load HINACT
    for i=1:length(H)
        H(i).values(H(i).values>0.99*HNOFLO)=NaN;
        H(i).values(H(i).values>0.99*HDRY  )=NaN;
    end
catch ME
    error('mfLab:maskH:NoData',...
        '%s can''t find HINACT.mat or can''t fine one of HDRY, HNOFLO',mfilename);
end


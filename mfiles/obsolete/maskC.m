function C=maskC(C)
% C = maskH(C)
% masks C(i).values using HINACT.mat saved by mf_setup.
%
% EXAMPLE:
% C = maskC(C);
% where C is the struct provided by readDat([basename '.HDS'];

% SEE ALSO:
%   maskBud readDat readMT3D readBud maskHC maskH
%
% TO 090315 100420 100718 100911 120409 120509 120524
%
% Copyright 2009-2012 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

if nargin<1, error('mflab:maskHC:insufficientInputArguments',...
    'maskC: 1 input argument required:\n');
end

try
    load CINACT
    for i=1:length(C)
        C(i).values(C(i).values==CINACT)=NaN;
    end
catch ME
    error('mfLab:maskC:NoMASK',...
        '%s can''t find CINACT.mat or can''t fine one of HDRY, HNOFLO',mfilename);
end

%% It seems stupid but MT3DMS does write extra records at some randdom STRESS PERIOD
% Probably because intertally there are small numerical differences so that
% the specified time is does not numericaly equal to the stres point end time.
% Anyway the difference is only about 1e-6.
% As a workaround, we remove these extra records here to prevent bothering
% the user.
C=C([1 find(diff([C.time])>1e-3)]);

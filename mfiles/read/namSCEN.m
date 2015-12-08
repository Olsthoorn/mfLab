function scenNr=namSCEN(PCKG,PCKGS,SCENS,varargin)
%NAMSCEN looks to see if PCKG is in PCKGS. If not, returns 0
%
% Example:
%    scenNr=namSCEN(PCKG,PCKGS,SCENS[,'exact'])
%
%    look to see of PCKG is in PCKGS if not return 0
%    if it is, return the correspondign SCEN Nr.
%
%    using in mf_setup to check if package PCKG is switched on in the NAM
%    worksheet.
%
% INPUTS:
%    PCKG is string
%    PCKGS a cell array with package names
%    SCEN an array with integers, same size as PCKGS
%
% TO 110807

i=strmatchi(PCKG,PCKGS,varargin{:});
if i==0,
    scenNr=0;
else
    scenNr=SCENS(i);
end
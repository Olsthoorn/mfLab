function newdem=mf_cleandem(dem)
%MF_CLEANDEM removes leading and trailing rows and columns with NaNs from dem
%
% Example:
%    newdem=mf_cleandem(dem)
%
%
% TODO: check, generalyze
%       used in mflab/examples/mf2007/Khetta_Jorf
%
% TO 110531

z=dem.z;

% find columns with NaNs
m=min(z);      Ie=find(~isnan(m)); Ie=Ie(1):Ie(end);

% find rows with NaNs
m=min(z,[],2); In=find(~isnan(m)); In=In(1):In(end);

newdem.eGr=dem.eGr([Ie Ie(end)+1]);
newdem.nGr=dem.nGr([In In(end)+1]);
newdem.em =dem.em( Ie);
newdem.nm =dem.nm( In);
newdem.z=dem.z(In,Ie);
newdem.BoundingBox=[[newdem.eGr(1); newdem.eGr(end)],[newdem.nGr(1); newdem.nGr(end)]];
newdem.Width =length(Ie);
newdem.Height=length(In);
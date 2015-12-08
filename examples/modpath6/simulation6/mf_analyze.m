%% Analyzing output of the model, modpath6 manual, SIMULATIION5, pp47ff
%% Backward multiple release endpoint simulation.
% TO 091011 091129 120413

%% make sure path is set to parent directory
d=dir('..');
if ~strmatchi('mf_analyzeALL.m',{d.name})
    error('%s: Can''t find file %s in the parent directory,\n%s\n',...
        mfilename,'mf_analyzeALL.m',fileparts(pwd));
else        
    theseLayers = 5;
    
    mf_analyzeALL;
end

%% Showing particles
pGrp = pGrp.getEndPoints('mp6.endp');
pGrp.plotStartPoints(1,'b.');
pGrp.plotEndPoints([],'r.');
view(3);

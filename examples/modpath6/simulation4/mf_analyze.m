%% Analyzing output of the model
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
pGrp = pGrp.getTsrPoints('mp6.tsr');
pGrp.plotTsrPoints();
view(3);

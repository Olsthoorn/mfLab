%% Analyzing output of MODFLOW/MODPATH SIMULATION2 from manual.
% TO 130221

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
%figure; hold on; set(gca,'xlim',gr.xGr([1 end]),'ylim',gr.yGr([end 1]));
pGrp = pGrp.getEndPoints('mp6.endp');

pGrp.dispEndPoints();

pGrp.endPointStatistics;
pGrp.plotEndPoints([],'k.');
pGrp.plotStartPoints([],'b.');

pGrp = pGrp.getPathLines('mp6.plines');
pGrp.plotPath('r');
view(3);

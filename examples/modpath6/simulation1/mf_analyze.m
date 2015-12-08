%% Analyzing output of MODFLOW/MODPATH SIMULATION1 from manual.
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
%pGrp = pGrp.getPathLines('mp6.plines');

%pGrp.dispEndPoints();

pGrp.endPointStatistics;
pGrp.plotStartPoints(2,'b.');
pGrp.plotStartPoints(3,'r.');
%pGrp.plotPath('r');
view(3);

% %% Showing particles
% figure; hold on; set(gca,'xlim',gr.xGr([1 end]),'ylim',gr.yGr([end 1]));
% pGrp = pGrp.getEndPoints('mp6.endp');
% pGrp.endPointStatistics;
% pGrp.plotEndPoints(1,'ko');
% pGrp.plotStartPoints(2,'b.');
% pGrp.plotStartPoints(3,'r.');
% 
% %%
% [endp,endpStat] = readEndPoints('mp6.endp');
% [tsr ,tsrStat ] = readTsrPoints('mp6.tsr');
% [path,pathStat] = readPath(     'mp6.plines');
% 
% %%
% figure; hold on; title('end points');
% plotStartPoints('mp6.endp','b.');
% plotEndPoints('mp6.endp','ro');
% plotTsrPoints('mp6.tsr','kp');
% plotPath('mp6.plines','g');
% view(3);
% 
% %%
% 
% figure; hold on; title('time series points');
% pGrp = pGrp.getTsrPoints('mp6.tsr');
% pGrp.tsrStatistics
% pGrp.plotTsrPoints();
% view(3);
% 
% %%
% figure; hold on; title('paths');
% pGrp = pGrp.getPathPoints('mp6.plines');
% pGrp.pathStatistics
% pGrp.plotPath('k');
% view(3);
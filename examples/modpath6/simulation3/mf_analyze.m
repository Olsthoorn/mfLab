%% Analyzing output of MODFLOW/MODPATH SIMULATION4 from manual.
% TO 130221

%% make sure path is set to parent directory
d=dir('..');
if ~strmatchi('mf_analyzeALL.m',{d.name})
    error('%s: Can''t find file %s in the parent directory,\n%s\n',...
        mfilename,'mf_analyzeALL.m',fileparts(pwd));
else     
    theseLayers = 5;  % contour layer 5 in mf_analyzeALL.m
    
    mf_analyzeALL;
end

%% Showing particles

%figure; hold on; set(gca,'xlim',gr.xGr([1 end]),'ylim',gr.yGr([end 1]));
% Then on the same axis plot the path lines and the their starting points

pGrp = pGrp.getEndPoints('mp6.endp');

%pGrp.dispEndPoints();

pGrp.endPointStatistics;
%pGrp.plotEndPoints(1,'ko');
pGrp.plotStartPoints(2,'b.');
pGrp.plotStartPoints(3,'r.');
pGrp = pGrp.getPathLines('mp6.plines');

pGrp.plotPath('b');
view(3);

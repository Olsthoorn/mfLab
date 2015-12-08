%% Analyzing model output

load('name'); load(basename); load('underneath');

%% Reading unformattted files
H=readDat([basename,'.hds']);  % just have a look (can also look in the LST file
C={}; Cs={};
C{1}= readMT3D('MT3D001.UCN');
Cs{1}=readMT3D('MT3D001S.UCN');

%% Get layer parameters used from workbook to construct a useful figure title

[LAYparnams,LAYparvals]=getLayers([basename '.xls'],gr);
aL     =LAYparvals(:,strmatchi('AL'  ,LAYparnams));
Rho    =LAYparvals(:,strmatchi('RhoB',LAYparnams));
Kd     =LAYparvals(:,strmatchi('SP1_1',LAYparnams));
Lambda =LAYparvals(:,strmatchi('RC1_1',LAYparnams));
R=1+Rho.*Kd./peff;

%% Then plot the results for the 4 layers in separte subplots
[~,~,~,~,NCOMP]=getExcelData(basename,'MT3D','V','NCOMP');

for iComp=1:NCOMP
    figure;
    
    it=length(C{iComp});
    for iLay=1:gr.Nlay
        subplot(gr.Nlay,1,iLay,'nextplot','add','xlim',gr.xGr([1 end]),'ylim',[0 C0],'xgrid','on','ygrid','on');
        
        xlabel('xGr [m]');
        ylabel('Relative conc');
        
        s=sprintf('Mt3D %s, species %d, layer=%d, time=%.1f d, aL=%g m, R=%g, Lambda=%g/d ',...
            basename,iComp,iLay,C{iComp}(end).time,aL(iLay),R(iLay),Lambda(iLay));
        
        title(s);
        plot(gr.xc,C{iComp}(it).values(1,:,iLay),[mf_color(iLay) '-']);

    end
end

% That's all, compare with the manual
% TO 091202
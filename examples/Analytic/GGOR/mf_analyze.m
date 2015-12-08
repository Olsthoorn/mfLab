%% Analyzing the simulation reusult of the shallow Dtuch top system.
% TO 100823
%
% Copyright 2009 2010 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

clear variables; close all;

load name
load(basename)
load underneath

%verification;
 
%% Build the model grid

[xGr,yGr,Z,xm,ym,Zm,DX,DY,DZ,Nx,Ny,Nz]=modelsize3(xGr,yGr,zGr);

%% Get heads of first layers 

H=readdat([basename '.HDS']);

dDYDX=(IBOUND(:,:,1)~=0 & IBOUND(:,:,1)~=iDITCH1).*(Dy*Dx);  % Area of active cells not being ditch
dDYDX=dDYDX./(sum(dDYDX,2)*ones(size(Dx)));                  % total area of cross section

h1=zeros(length(P),length(H));
h2=zeros(length(P),length(H));
for i=1:length(H)
    h1(:,i)=sum(H(i).values(:,:,1).*dDYDX,2);
    h2(:,i)=sum(H(i).values(:,:,2).*dDYDX,2);
end

%% test


try
    axes(ax(1));
    h=findobj(get(gcf,'children'),'Tag','legend');
    s=get(h,'string');    
catch ME
    s={};
    figure; hold on; xlabel('time [d]'); ylabel('head [m]');
    title(sprintf('Simulation of dynamic groundwater heads in %d parcels in %s',length(P),basename));
end

t=tne(:,1)'; 
plot(t,h1( 1,:),'r' , 'linewidth',0.75); s=[s, 'MFLOW h(1)'];
plot(t,h2( 1,:),'g.' , 'linewidth',0.75); s=[s, 'MFLOW h(2)'];

datetick('x');
legend(s);

%%
% axes(ax(2));
% 
% h=findobj(get(gcf,'children'),'Tag','legend');
% s=get(h,'string');
% 
% t=tne(:,1)'; 
% plot(P.b-xm(end:-1:2),H1(end).values(:,end:-1:2),'r.' , 'linewidth',0.75); s=[s, 'MFLOW h(1)'];
% plot(P.b-xm(end:-1:2),H2(end).values(:,end:-1:2),'r+' , 'linewidth',0.75); s=[s, 'MFLOW h(2)'];
% 
% figure; hold on;
% plot(xm,H1(end).values,'r.' , 'linewidth',0.75); s=[s, 'MFLOW h(1)'];
% plot(xm,H2(end).values,'r+' , 'linewidth',0.75); s=[s, 'MFLOW h(2)'];
% 
% 
% legend(s);

%% === water budgets ===== water budgets ===== water budgets =====

    
[glg,gvg,ghg]=getGXG(h1,tne(:,1),2);

for iP=1:length(P);
    P(iP).GHG=ghg(iP);
    P(iP).GVG=gvg(iP);
    P(iP).GLG=glg(iP);
end

if size(h1,1)>1   
    %% plot GxG from model and analytic
    leg=[];
    figure; hold on; grid on;
    plot([P.GLGDBF],'r--'); leg{end+1}='GLGDBF';
    plot([P.GVGDBF],'g--'); leg{end+1}='GVGDBF';
    plot([P.GHGDBF],'b--'); leg{end+1}='GHGDBF';
    plot([P.GLG],'r');      leg{end+1}='GLG';
    plot([P.GVG],'g');      leg{end+1}='GVG';
    plot([P.GHG],'b');      leg{end+1}='GHG';
    xlabel('Cross section number');
    ylabel('head (above datum)');
    grid on
    title('GxG of the computed cross sections (GGOR tool, MODFLOW)');
    legend(leg);
end

return

fprintf('press RETURN to continue\n');
pause;

%% Compute budgets

userLabels={'STORAGE' 'WELLS' 'DRAINS' 'HEADDEPBOUNDS' 'RECHARGE'}; % needed for water balance

compact=1;

if compact>0
    B=readBud1(['-' basename '.bgt'],userLabels); % selected cross sections for water balance
else
    B=readBud1([    basename '.bgt'],userLabels); % selected cross sections for water balance
end

                           %% Narrative water budgets
% The budget file has the flow terms for each flow proejct (see labels) in
% 4D, that is for each time we have a 3D array of values for each flow
% process. This is an amazing amount of data.
% First we will sum the data flow terms for every cross sections, so that
% we end up with the total of DRN, GHB ect for every cross section and time
% period. We can produce a water balance for every individual cross section
% over time from this.
% We will then be able to combine cross sections arbitrarily if we so
% desire. We may for instance compute the overall running water balance for
% the entire area by summing over all cross sections weighted by their
% area. We could do so for every process separately, to find the total
% runoff (= total DRN), total recharge total ditch discharge (= GHB) and
% total seepage (= given WEL). The combination of arabitrary subsaries into
% different subareal totals can be done, if these subareas have indices in
% the database that can be used to recognize the subareas. The data base is
% contained in the structmatrix P.

%% Terms for the water budget
% The terms of interest in this model are
% RCH  -- recharge
% GHB  -- flow from the ditches
% DRN  -- overflow (surface runoff)
% WEL  -- fixed seepage flow injected in second aquifer
% STO  == storage during time step
% to compare these and check the water balans, we may compute the values
% per unit of ground surface m/s or mm/d

%% To get the discharges of the individual cross sections per unit area
%  in m/d we first sum all flows over the entire cross section and then
%  have to divide the values by the surface area of the cross section

NSec=length(B(1).rows);    % number of cross sections in budget file
                           % same as length ISec (selected cross sections)

%% Get active area of the cross sections from IBOUND<>0 in layer 1

a=DY.*(xGr(sum(IBOUND(:,:,1)~=0,2)+1)'-xGr(1));

%% Allocate memory to store budget terms.
% NSec = cross section, Nt is time

Nt=size(tne,1);
QRCH=zeros(NSec,Nt);   % flow form recharge
QGHB=zeros(NSec,Nt);   % flow through the general head boundaries (to and form the ditches)
QDRN=zeros(NSec,Nt);   % flow through the drains at ground surface
QWEL=zeros(NSec,Nt);   % prescribe flow for the seeapage from or to the regional aquifer
QSTO=zeros(NSec,Nt);   % storage change during the time step

% Fill the arrays using the flow terms from the budget file.

idrn = strmatchi('DRAINS'       ,B(1).label);
irch = strmatchi('RECHARGE'     ,B(1).label);
ighb = strmatchi('HEADDEPBOUNDS',B(1).label);
isto = strmatchi('STORAGE'      ,B(1).label);
iwel = strmatchi('WELLS'        ,B(1).label);

% Sum over columns and layers to get totals for all cross sections as a
% funcion of time

if compact>0   % compact form of budget arrays X=[I V], [globindex value] 
    dims=[length(B(1).rows) length(B(1).cols) length(B(1).lays)];

    for i=1:length(B)
        QDRN(:,i)=sum(sum(mf_expand(B(i).term{idrn},dims),2),3); % catch end
        QRCH(:,i)=sum(sum(mf_expand(B(i).term{irch},dims),2),3); % catch end
        QGHB(:,i)=sum(sum(mf_expand(B(i).term{ighb},dims),2),3); % catch end
        QSTO(:,i)=sum(sum(mf_expand(B(i).term{isto},dims),2),3); % catch end
        QWEL(:,i)=sum(sum(mf_expand(B(i).term{iwel},dims),2),3); % catch end
    end    
else
    for i=1:length(B)
        QDRN(:,i)=sum(sum(B(i).term{idrn},2),3); % catch end
        QRCH(:,i)=sum(sum(B(i).term{irch},2),3); % catch end
        QGHB(:,i)=sum(sum(B(i).term{ighb},2),3); % catch end
        QSTO(:,i)=sum(sum(B(i).term{isto},2),3); % catch end
        QWEL(:,i)=sum(sum(B(i).term{iwel},2),3); % catch end
    end
end

% compute these discharge in m/d by dividiing with the cross section area
QDRN=QDRN./(a*ones(1,Nt));
QRCH=QRCH./(a*ones(1,Nt));
QGHB=QGHB./(a*ones(1,Nt));
QSTO=QSTO./(a*ones(1,Nt));
QWEL=QWEL./(a*ones(1,Nt));

%% Total these flows to see of the water balance is zero (or almost)
QTOT=         QDRN +    QRCH +    QGHB +    QSTO +    QWEL;
QABS=0.5*(abs(QDRN)+abs(QRCH)+abs(QGHB)+abs(QSTO)+abs(QWEL));

%% Display these flows summed over all times for the cross sections
% to show that the water budget matches
% for all cross section of all times
fprintf('      QTOT      QWEL      QSTO      QDRN      QGHB      QRCH  [mm/d]');
display(1000*[...
    mean(QTOT,2),...
    mean(QWEL,2),...
    mean(QSTO,2),...
    mean(QDRN,2),...
    mean(QGHB,2),...
    mean(QRCH,2)]);

%% saving the waterbalance for all sections
labels=B(1).label;
timestamp=now;
legtxt={'UPWSPG','STO   ','DITCH ','RUNOFF','RCH   ','DWNSPG','STO   ','DITCH ','RUNOFF','EVTR  '};
labtxt={'QWEL>0','QSTO>0','QGHB>0','QDRN>0','QRCH>0','QWEL<0','QSTO<0','QGHB<0','QDRN<0','QRCH<0'};
save Balance labels legtxt labtxt QTOT QWEL QSTO QDRN QGHB QRCH P tne timestamp

%% The (relative) area of the selecred cross sections

A=[P.area]'; a = (A/sum(A))*ones(1,size(QTOT,2));

%% Show running water budget summed over all area-weigted cross sections 
%  selected from the database

POS=[max(0,sum(a.*QWEL,1));max(0,sum(a.*QSTO,1));...
     max(0,sum(a.*QGHB,1));max(0,sum(a.*QDRN,1));...
     max(0,sum(a.*QRCH,1))];
NEG=[min(0,sum(a.*QWEL,1));min(0,sum(a.*QSTO,1));...
     min(0,sum(a.*QGHB,1));min(0,sum(a.*QDRN,1));...
     min(0,sum(a.*QRCH,1))];
figure; hold on

% convert to mm/d while plotting
area(tne(:,1),POS'*1000);
area(tne(:,1),NEG'*1000);
plot(tne(:,1)',sum(a.*QTOT*1000,1),'w');

legend('UPWSPG','STO','DITCH','RUNOFF','RCH','DWNSPG','STO','DITCH','RUNOFF','EVTR');
xlabel(sprintf('%s - %s',datestr(tne(1,1),'yyyy'),datestr(tne(end,1),'yyyy')));
ylabel('mm/d');
title(sprintf('Running water budget % %s over %d cross sections representing %.0f ha in mm/d',...
    basename,NSec,sum([P.area])/1e4));
datetick('x',4); grid on;

%% Show the running balance of the top system

%dtopBalance

function [M,time] = massBudget(o,iComp,space)
% [M,time] = well.massBudget(iComp,space);   % compute massbudget of well over time
%    iComp = compoment nr (species nr)
%    space = [r dz], where r is the distance around the wells to be included
%       in the mass balance and dz the space above and below the screens. Total
%       space for which mass balance is computed is +/-r around well and +/- z
%       above and below the screen.
%       STCONC= starting concentrations cell array whre STCONC{iComp} is the array for compoment iComp
%
%   The change of mass is Qin cin dt or -Q cout dt where cin is the injected concentration
%   and cout the extracted concentration. This equals the change of mass in
%   the model over the same time increment. Starting with the STCONC at t=0
%   and the computed concentrations for all other times, we can compute the
%   total mass change between different times directly from the computed
%   concentrations in the model. We compute the total mass for a range
%   around the well specified by space. If space is large it may thus
%   specify the entire model. The resuls are correct if the concentrations
%   at the model boundaries remain unchanged. This is the case if ICBUND is
%   negative (fixed concentrations), or at least when there is no
%   concentration gradient at the boundary. In all other cases this
%   boundary mass flux must be taken into account, which is not done by
%   this procedure.
%
%
%  TO 121122 (Well Limburg)

load('name'); load(basename);

B = readBud([basename '.BGT']);
[LAYnams,LAYvals] = getLayers(basename,B(1).NLAY);

% reaction package on?
[Nnum,Ntxt]=xlsread(basename,'NAM','','basic');
isonRCT = Nnum(strmatchi('RCT',Ntxt(:,1)),3)>0;

% which type of sorption?
[MT3Dnams,MT3Dvals] = getExcelData(basename,'MT3D','Vertical');
ISOTHM = MT3Dvals(strmatchi('ISOTHM',MT3Dnams),1);

ISOTHM = ISOTHM *  isonRCT;

rhob = LAYvals(:,strmatchi('RHOB',LAYnams));
SP1   = gr.const(LAYvals(:,strmatchi(sprintf('SP%d_1',iComp),LAYnams)));
SP2   = gr.const(LAYvals(:,strmatchi(sprintf('SP%d_2',iComp),LAYnams)));

RHOB = gr.const(rhob);

if ~numel(space) == 2, error('%s: arg 5 must be a two number vector [r dz]'); end

if ~iscell(STCONC), STCONC = {STCONC}; end %#ok

%% Get the dissolved and sorbed concentrations
C  = readMT3D(sprintf('MT3D%03d.UCN',iComp));
% if exist(sprintf('MT3D%03dS.UCN',iComp),'file')
%     CS = readMT3D(sprintf('MT3D%03dS.UCN',iComp));  % sorbed conc
% else
%     CS = C;
%     for it = 1:length(C), CS(it).values(:) = 0; end
% end

%% Set Cout of wells
o = o.setCout(C,iComp,B);

%% for all wells having negative input concentrations (denoting index of the
% recirculation well, i.e. the cell whose concentration is used,
%  get the concentration of the output.
for iw = 1:numel(o)
    for it=1:length(o(iw).Q)
        if o(iw).C(iComp,it)<0,
            o(iw).C(iComp,it) = C(it).values(abs(o(iw).C(iComp,it)));
        end
    end
end

%% (rectangular) Zone around each well using its well number as zone number
r = space(1); dz=space(2);
ZONE = zeros(gr.size);
for iw = 1:numel(o)
    z1 = max(o(iw).z);
    z2 = min(o(iw).z);
    ZONE(gr.YM>o(iw).y(1)-r & gr.YM<o(iw).y(1)+r &...
         gr.XM>o(iw).x(1)-r & gr.XM<o(iw).x(1)+r &...
         gr.ZMlay>z2-dz     & gr.ZMlay<z1+dz) ...
     = o(iw).nr;
end

time = cumsum(o(1).Dt);
NT   = length(o(1).Q);

% Total mass in area surrounding the well
% since we don't know the sorbed conc at t=0, we compute everything
% relative to time 2

nCols = 5;

M = NaN(NT,4*iw);
for iw=numel(o):-1:1
    IDX = find(ZONE==o(iw).nr);
    
    for it = NT:-1:1
        % Injected and extracted mass
        if o(iw).Q(it)>0
            M(it,nCols*(iw-1)+1) = o(iw).Q(it).*o(iw).C(   iComp,it) * o(iw).Dt(it);
        else
            % Extracted mass
            M(it,nCols*(iw-1)+2) = o(iw).Q(it).*o(iw).Cout(iComp,it) * o(iw).Dt(it);
        end
        
        % Self computed dissolved mass
        M( it,nCols*(iw-1)+3) = sum(PEFF(IDX).*gr.Vlay(IDX).*C(it).values(IDX));
        
        % Self computed sorbed mass
        M( it,nCols*(iw-1)+4) = sorbedMass(C(it).values(IDX),ISOTHM,gr.Vlay(IDX),RHOB(IDX),SP1(IDX),SP2(IDX));
        
        % MT3DMS computed sorbed mass
        %M( it,nCols*(iw-1)+5) = sum(RHOB(IDX).*gr.Vlay(IDX).*CS(it).values(IDX));
        
        % Dissolved + Sorbed
        M( it,nCols*(iw-1)+5) = M(it,nCols*(iw-1)+3)+M(it,nCols*(iw-1)+4);
    end
    M(:,nCols*(iw-1)+1) = cumsum(M(:,nCols*(iw-1)+1));
    M(:,nCols*(iw-1)+2) = cumsum(M(:,nCols*(iw-1)+2));

    % Self computed initial mass 
    % MSI = sorbedMass(STCONC{iComp}(IDX),ISOTHM,gr.Vlay(IDX),RHOB(IDX),SP1(IDX),SP2(IDX));
     %M(:,nCols*iw-0) = [ MSI; diff(M(:,nCols*iw-0,1,1))];
     %M(:,nCols*iw-1) = [ MSI; diff(M(:,nCols*iw-1,1,1))];
end

%% Print zone budgets for all wells

fprintf('1) time in user units at end of each time step or stress period.\n');
fprintf('2) For each well the following columns are printed:\n');
fprintf('     QDT when the well is injecting\n');
fprintf('     QDT when the well is extracting\n');
fprintf('     change of Mass around well in area +/-r and +/- Dz\n');
fprintf('A = Injected [kg]\n');
fprintf('B = Extracted  [kg]\n');
fprintf('C = Dissolved [kg]\n');
fprintf('D = Sorbed    [kg]\n');
fprintf('E = Sorbed (MT3DMS/SEAWAT)\n');
fprintf('%12s','time');
hdr = repmat({o(iw).name},[1,nCols]);
for iw = 1:numel(o),
    fprintf(' %10s_A %10s_B %10s_C %10s_D %10s_E',hdr{:});
end
fprintf('\n');
fmt = ['%12g' repmat(repmat(' %12g',[1,nCols]),[1,numel(o)]) '\n'];

fprintf(fmt,[time' M]');

fprintf('\n\n');


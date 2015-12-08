function [h,Qfix,Qnod,QstoSS,QCau,Qz,z,horizon,cropNm,SS,Theta,Kh,Acrop,Intercepted]=unsat1(TPE,dataWorkBook,varargin)
%UNSAT1t a 1D block-centred transient FD model for the unsaturated zone.
% One column is simulated. Seepage is added to the lower cells,
% precipitatio to the top cells higher than the defined infDepth. The
% evapotranspiration is subtracted from the root zone as defined by the
% root-zone depth.
% The input P an E are filtered by interception before applying the so
% reduced P and E to the model.
% There is no fixed head defined, but all cells below the defined ditch
% elevation are coupled to the ditch through their portion of the drainage
% resistance.
% The initial moisture is obtained from the initial suction head. The
% initial suction head is assumed to be equal to the elevation above the
% ditch.
% An arbitrary number of horizons may be specified in the params sheet by
% using the soil name as label and the bottom elevation of the horizon as
% value. The deepest horizon will automatically be extended to the bottom
% of the profile.
% Notice that z is positive, in m, the depth below ground surface.
% The suction head will also be in m, it is easily converted to pF by
% pF = log10(-h*100);
% Notice that suction head is positive where p<atmospheric, if
% p>atmospheric we talk about pressure head.
%
% Example:
%    [h,Qt,Qz,Qs]=unsat1(TPE,dataWorkbook);
%
% INPUT:
%  TPE is the array with three columns, [time P E], P and E in m/d
%  z   = z-coordinate of mesh/grid; a vector
%
% The Excel data workbook stores the parameters used in sheet 'params'
% It stores the soil data in sheet 'soilData'.
% The soilData stem from Dingman (2002), table 6-1, page 233.
% It stores the parameters for 11 soils published by Clap and Hornberger (1974)
% for the Campbell analytical soil relations. See Dingman (2002) eq. 6-12 and 6-13.
%
% Use (a copy of) the Excel workbook ClapHornberger78SoilData.xls to adapt
% the data to your specific case.
%
% Notice that all parameters are in m and days not in cm and days
%
%
% OUTPUT
%  h(Nz,1)    suction head (cm) is is -p/(rho g)
%, Qt(Nz,Nt)    computed total cell balance during timestep it (m/d)
%  Qz(Nz-1,Nt)  hor.  cell face flow in z-direction downward positive(m/d)
%  Qs(Nz,Nt)    storage change of node during timestep it (m/d)
%
% See also: fmd2t fdm2c fdm2ct fdm3 fdm3t
%
% TO 150911 991017  TO 000530 001026 070414 080301

% Copyright 2015 Theo Olsthoorn, without any warranty
% under free software foundation GNU license version 3 or later

if numel(varargin)>0
    sheetNm = varargin{1};
else
    sheetNm = 'params';
end
 
    %% Variable simulation input
    t  = TPE(:,1);
    dt = diff(t); dt=[dt(1); dt]; % to compute also the end of the first day !!
    Nt = numel(dt);
    P0 = TPE(:,2);
    E0=TPE(:,end);
    
    %% Get simulation parameters and soil property data
    
    % Default data workbook
    if nargin<3, dataWorkBook = 'ClapHornberger78SoilData'; end
    
    % Parameters and sol data
    [dataLbls,data] = getExcelData(dataWorkBook,sheetNm,'V');

    % shortcut for extracting paramters from worksheet
    dat = @(v) data(strmatchi(v,dataLbls,'exact'),1);
    
    %% Model configuration and grid
    hGround   = dat('maaiveld');  % [m] NAP ground surface 
    dCover    = dat('dCover');   % [m]   working thickness of cover layer
    dCell     = dat('cellThickness');   % [m] dikte modelcellen
    
    iPond     = 1; % cell number of pond on ground surface
    
    z  = hGround- (0:dCell:dCover); z =z(:);
    dz = abs(diff(z));
    zm = 0.5*(z(1:end-1)+z(2:end));
    Nz = numel(dz);

    %% Divison of infiltration and evapotranspiration over cells
%    infDepth  = dat('infDepth');  % [m] bottom of infiltration zone, top at z=0
    rootDepth = dat('rootDepth'); % [m] bottom of root zone, top = at z=0
    rootZone= zm>z(1)-rootDepth;
    fInf  = zeros(Nz,1);  fInf(iPond)=1;
    fRoot = zeros(Nz,1);
    fRoot(rootZone) = dz(rootZone)./sum(dz(rootZone)); % sums to 1.0

    %% parcel and ditch
    hDitch    = dat('hDitch');    % [m] elevation of water level in ditch    % model setup
    Lin       = dat('Lin');         % [m] width of field or parcel
    Lex       = dat('Lex');       % [m] width of subfield drainage system
    wIn       = dat('wIn');       % [d] entry resistance of ditches    
    wEx       = dat('wEx');       % [d] exit resistance of ditches
    omega     = dat('omega');     % [m] wet circumference of ditches
    kCover    = dat('kCover');    % [m/d]  conductivity of cover layer
    fixedH    = dat('fixedH');    % [ m ]  fixed pressure head at bottom of model or NaN
    
    %% Drainage
    gammaRunoff = dat('gammaRunoff'); % drainage resistance surface runoff
    gammaDrain  = dat('gammaDrain');  % drain resistance use inf if no drain
    dDrain      = dat('dDrain');      % drain depth
    
    if dDrain<0, error('dDdrain must be > 0 !, not %g m',dDrain); end

    sigmaPond = dat('sigmaPond'); % [ m ]  std of ground surface unevenness, for ponding
    Cpond = @(h) 0.5*erfc(-h/sigmaPond);                  % [-] storage coefficient ponding
    Wpond = @(h) 0.5 * sigmaPond * ierfc(-h/sigmaPond,1); % [m] storage in ponding
    
    %% Soil data choice
    Staring   = dat('Staringreeks')     ~=0; % use Staring reeks or USA soils with Van Genuchten relations
    usaCorey  = dat('Hornberger Series')~=0; %use Hornberker parameters and Corey formulas + USA soils
    
    if Staring, usaCorey=false; end

    %% Computional data
    maxOuter  = dat('maxOuter');  % [m] maximum number of outer iterations
    hClose    = dat('hClose');    % [m] head closure criterion
    upWind    = dat('Upwind')~=0; % use upwind FD scheme for better advection
    
    implicity = dat('implicit');  % [-] implicitness    
    % verify implicitness
    if implicity<0.5 || implicity>1
        error('implicitness must be >=0.5 and <=1, not %g',implicity);
    end
        
    %% Interception and seepage
    maxIc     = dat('maxIc');     % [m] maximum interception
    seepage   = dat('seepage');   % [m/d] upward seepage (downward if <0)    
    
    % Filter recharge and potential evapotranspiration
    Intercepted = max(0,min(min(maxIc,P0),E0));
    P = P0-Intercepted;
    E = E0-Intercepted;
    Pmax = max(0,P);

    %% Soil horizons and soil properties
    if usaCorey
        [soilLbls,soilData,~,soilNames] = getExcelData(dataWorkBook,'usaCorey','H');
    else
        [soilLbls,soilData,~,soilNames] = getExcelData(dataWorkBook,'Staringreeks','H');
    end
        
    I = find(ismember(dataLbls,soilNames));
    horizon = dataLbls(I);
    for i=numel(I):-1:1
        depth(i)  = data(I(i),1);
        if depth(i)==0
            horizon(i)=[];
            depth(i)  =[];
        end
    end
    if numel(horizon)==0
        error(['No soil horizons are specified in the unsat sheet of file %s\n', ...
               'Specify at least one soil as label and one depth in sheet unsat\n',...
               'See recognized soils in sheet soils'],dataFile);
    else % sort soils and depths
        if any (depth<=0), error('horizon depths must all be > 0, except for Pond10cm'); end
        [depth,I] = sort(depth,'descend');
        horizon   = horizon(I);
        depth(1)  = z(1)-z(end);    % largest z covers rest of profile by default
    end
       
    %% Crop, regulates root-water uptake reduction/stress
    [cropLbls,cropData,~,cropNames] = getExcelData(dataWorkBook,'cropData','H'); 

    I  = find(ismember(dataLbls,cropNames));
    cropNm = dataLbls(I(data(I,1)>0));

    cropPar = @(cropNm,par) cropData(strmatchi(cropNm,cropNames,'exact'),strmatchi(par,cropLbls,'exact'));

    if numel(cropNm)~=1
        error('Specify one and only one a crop name in sheet cropData of workbook %s\n',dataWorkBook);
    end
    
    % Only uses Feddes drought stress, not saturation stress
    h3 = cropPar(cropNm{end},'h3');  % [ cm ]alternative h3Low and h3High, set in crop sheet
    h4 = cropPar(cropNm{end},'h4');  % [ cm ];
    % Units conversion
    h3 = h3/100; % [ m ];
    h4 = h4/100; % [ m ];
    
    %% initialize soil parameters
    if usaCorey
        %% usa soils and Corey formulas are used see Dingman (2001, table 1)
        KsCor = zeros(Nz,1);
        peff = zeros(Nz,1);
        bCor = zeros(Nz,1);
        psiAE= zeros(Nz,1);

        % shortcut to extract soi parameter for given soil
        soilPar = @(soilNm,par) soilData(strmatchi(soilNm,soilNames,'exact'),strmatchi(par,soilLbls,'exact'));

        for is=1:numel(horizon)
            KsCor(zm>z(1)-depth(is)) = soilPar(horizon{is},'KsCor'); % [cm/d]
            peff( zm>z(1)-depth(is)) = soilPar(horizon{is},'peff');  % [ - ]
            bCor( zm>z(1)-depth(is)) = soilPar(horizon{is},'bCor');  % [ - ]
            psiAE(zm>z(1)-depth(is)) = soilPar(horizon{is},'psiAE'); % [ cm ]
            % Units conversion
            KsCor = KsCor/100; % [ m/d ]
            psiAE = psiAE/100; % [  m  ]
        end

        %% Parameters we need to make theta-psi and Ss-psi continuous
        % used in theta-psi relation and Ss relation
        delta        = 0.25;  % [-] small value to switch to assympothic curve
        psiStar      = psiAE * (1+delta); % psiStar [m]
        thetaStar    = peff  .* (1+delta) .^ (-1./bCor);
        dthDpsiStar  = peff./bCor ./ psiAE .* (psiStar./psiAE).^(-1./bCor-1);
        lambda       = (peff-thetaStar)./dthDpsiStar;  % in cm
         
    else
        %% Staringreeks + Van Genuchten parameters are used
        KsvG   = zeros(Nz,1);
        thetar = zeros(Nz,1);
        thetas = zeros(Nz,1);
        LvG    = zeros(Nz,1);
        nvG    = zeros(Nz,1);
        alpha  = zeros(Nz,1);

        % shortcut to extract soil parameter for given soil
        soilPar = @(soilNm,par) soilData(strmatchi(soilNm,soilNames,'exact'),strmatchi(par,soilLbls,'exact'));

        for is=1:numel(horizon)
            KsvG(  zm>z(1)-depth(is)) = soilPar(horizon{is},'KsvG');   % [ cm/d]
            thetar(zm>z(1)-depth(is)) = soilPar(horizon{is},'thetar'); % [ - ]
            thetas(zm>z(1)-depth(is)) = soilPar(horizon{is},'thetas'); % [ - ]
            LvG(   zm>z(1)-depth(is)) = soilPar(horizon{is},'LvG');    % [ - ]
            nvG(   zm>z(1)-depth(is)) = soilPar(horizon{is},'nvG');    % [ - ]
            alpha( zm>z(1)-depth(is)) = soilPar(horizon{is},'alpha');  % [1/cm]

            % Units conversion
            KsvG = KsvG/100;    % [m/d]
            alpha= alpha * 100; % [ 1/m] 
        end
        mvG = 1 - 1./nvG;
    end
    
    
    %% ========= With all parameters set, build the model ==================
    
    %% Initialize or allocate
    Cz   = NaN(Nz-1,1);
    Cs   = NaN(Nz,1);
    Cp   = NaN(Nz,1);
    K05  = NaN(Nz-1,1);
    Ccau = NaN(Nz,1);
    W    = zeros(Nz,1);

    % Save for later interpretations
    FQ      = zeros(Nz,1);   % Fixed inflow during time step = injection
    
    % For output and water budget
    h       = NaN(Nz, Nt); % at end of time step
    Theta   = NaN(Nz, Nt); % at end of time step
    Acrop   = NaN(Nz, Nt); % during time step
    Kh      = NaN(Nz, Nt); % of cell
    SS      = NaN(Nz, Nt); % during time step
    Qfix    = NaN(Nz, Nt);  % given inflow during time step
    Qnod    = NaN(Nz, Nt);  % net input of node during time step A*h
    QstoSS    = NaN(Nz, Nt); % water into storage during time step computed using Ss
    QstoTheta = NaN(Nz, Nt); % water into storage during time step computed using theta
    QCau    = NaN(Nz, Nt);  % Cauchy inflow during time step
    Qz      = NaN(Nz-1,Nt);  % Qz during time step
    K5      = NaN(Nz-1, Nt); % between cells
    Pond    = NaN(1,Nt);  % water height in ponds on ground surface

        %% node numbering and neighbor numbering
    Nodes = (1:Nz)';               % Node numbering
    Ibot  = Nodes(2:end,:);     % Top    neighbor node numbers
    Itop  = Nodes(1:end-1,:);   % Bottom neighbor node numbers
    %% Loop through time
    fprintf('Outer iterations for each time step, hClose=%gm:\n',hClose);
   
    hCauchy    = hDitch-zm;   % hydrostatic head boundary
    %hCauchy(zm>=hDitch) = 0;  % seepage face boundary equals p=0 or h=0
    strtH      = hCauchy;
    ht         = strtH;
    htEnd      = strtH;
    htOld      = strtH;
    
    iDrain = find(zm <= z(1) -dDrain,1,'first');

    Iact    = true(Nz,1);
    
    if ~isnan(fixedH), Iact(end)=false; strtH(end)=fixedH; end %#ok
        
    iCount=0; mxCount=10;
    for it=1:Nt
        fprintf('% d',it);
        dtau = dt(it);
        tau  = t(it)-dt(it);
        strtH = htEnd;
        
        FQ(:)   = P(it)*fInf - E(it)*fRoot.*aCrop(ht); % moet nog verdampingsredutie in  
        FQ(end) = seepage;

        % Outer loop because of non-linearity of system
        htOld(:)   = strtH;
        
        while t(it)-tau > eps;
            for iOuter = 1:maxOuter
                % resistances and conducctances
%                ht(1:2)     = RungeKutta(ht(1:2),dt(it));
%                htEnd(1)    = strtH(1) + (ht(1)-strtH(1))/implicity;
%                W(1)        = (Wpond(htEnd(1)) - Wpond(strtH(1)))/dt(it);
%                W(1)        = dz(1) * (h2theta(htEnd(1)) - h2theta(strtH(1))) / dt(it) + W(1);

                [Cz,K05,k] = Cond_z( ht);
                Ccau     = Ccauchy(ht);

                %% System matrix and main diagonal                
                A    =  sparse([Itop;Ibot], [Ibot;Itop], -[Cz;Cz], Nz,Nz,3*Nz);
                Adiag= -sum(A,2); 

                Cs(:)     = dz .* Ss(ht) / implicity;
                Cp(:)     = Cpond(ht(1));

                fprintf('dtau = ');
                while true;
                    fprintf(' %g',dtau);
                    RHS       = FQ + Ccau.*hCauchy + (Cs+Cp).*strtH/dtau + [0; K05] -[K05; 0];
                    A1        = spdiags(Adiag + Ccau + (Cs+Cp)/dtau,0,A);
                    ht(Iact)  = A1(Iact,Iact)\ (RHS(Iact) - A1(Iact,~Iact) * strtH(~Iact));
                    if ht(1)>z(1) + 2*Pmax  % ponds should not be deeper than reasonable
                        dtau=dtau/2;                    
                    else
                        break;
                    end
                end
                fprintf('\n');

                htEnd(:)  = strtH + (ht-strtH)/implicity;

                % proceed to end of time step using implicitness implicity
                if iOuter>1
                    dh    = sqrt(mean((htEnd-htOld).^2));
                    htOld  = htEnd;
                    iCount=iCount+1;
                    %fprintf('% g',dh);

                    % check and finish if no change of head between outer iterations
                    if dh<hClose || iOuter == maxOuter
                        break
                    end               
                end
            end
            tau = tau + dtau; % update current time
            dtau=dtau * 2;    % try lager time step
        end
        
        %% Compute flows only when final heads are known
        h(:,    it) = htEnd;
        Theta(:,it) = h2theta(htEnd);
        Acrop(:,it) = aCrop(h(:,it));
        Kh(   :,it) = k;
        K5(   :,it) = K05;
        Pond( it)   = W(1);

        Qz(:,it) = -Cz.*diff(ht)-K05;

        % Net outflow of the node
        Qnod(:,it) = spdiags(Adiag,0,A)*ht- [0; K05] + [K05; 0];

        % Fixed (given) inflow
        Qfix(:,it) = FQ;

        % Storage into the node
        QstoSS(   :,it) = Cs.*(ht-hCauchy)/dt(it);  % Cs includes implicity, so use ht not htEnd
        QstoTheta(:,it) = dz.*(Theta(:,it) - h2theta(strtH))/dt(it);

        % Net inflow through boundary Cauchy boundary (drainage
        % resistance)
        QCau(:,it)= Ccau .*(hCauchy - ht);
        fprintf('\n');

    end
    
    fprintf('Done, saving all data from unsat1 to unsat1.mat\n');
    
    % save results for use outside the model
    wie = who;
     save('unsat1',wie{:});
%      Acrop E Intercepted Kh KsCor KsvG L LDr LvG P ...
%          QCau Qfix Qnod Qsto Qz ...
%          S Staring TPE Theta alpha bCor cropNm ...
%          dCell dCover dataWorkBook ...
%          dthDpsiStar dz fInf fRoot gammaEx gammaIn h h3 h4 hCauchy ...
%          hClose hDitch  hGround horizon implicity ...
%          infDepth infZone kCover kD lambda maxIc ...
%          mvG nvG omega peff psiAE psiStar rootDepth rootZone satZone ...
%          seepage thetaStar thetar thetas usaCorey wEx wIn z zm 

    function theta = h2theta(h)        
        % compute theta from given pressure head 
        Nh = numel(h);
        if usaCorey
            theta   = peff(1:Nh);
            I = find(h < -psiStar(Nh)); 
            theta( I) = peff( I) .* (-h(I)./psiAE(I)).^(-1./bCor(I));
            theta(~I) = peff(~I) - (peff(~I) - thetaStar(~I)) .* exp((-h(~I)-psiStar(~I))./lambda(~I));        
        else % Staring
            theta = thetas(1:Nh);
            I = find(h<0);
            theta(I) = thetar(I) + (thetas(I) - thetar(I)) ./ ((1+(alpha(I).* -h(I)).^nvG(I)).^mvG(I));
        end
    end

    function ss = Ss(h)
        % requires a vector of length dz
        ssElastic = 0 * 1e-5;  % [1/m] ! 
        delta     = 0.02;  % characteristic ponding depth, exponential
        Nh = numel(h);
        ss = ssElastic * ones(Nh,1);
        if usaCorey
            I = find(h < -psiAE); % only if unsaturated
            ss( I)      = ss(I) + peff(I)./bCor(I)./psiAE(I).*(-h(I)./psiAE(I)).^(-1./bCor(I)-1);
            if any(~I)
                ss(~I)  = ss(~I) + dthDpsiStar(~I) .* exp((-h(~I)-psiStar(~I))./lambda(~I));
            end
        else % Staringreeks, + Van Genuchten
           I = find(h < 0); % only if unsaturated
           ss(I) = ss(I) + mvG(I).*nvG(I).*alpha(I).*(thetas(I)-thetar(I)).* ...
                (1+(alpha(I).*(-h(I))).^nvG(I)).^(-mvG(I)-1).*(alpha(I).*(-h(I))).^(nvG(I)-1);
        end
        if Nh>1, SS(:,it) = ss; end % save storage coefficient for later interpretation       
    end

    function k = K(theta)
        I = 1:numel(theta);
        if usaCorey
            % compute unsaturated hydraulic conductivity using Dingman (2002,
            % equation 6-13b with 6-14
            S = min(1,theta./peff(I));
            k = KsCor(I) .* S.^(2*bCor(I)+3);
        else % Staringrreeks + Van Genuchten
            S  = min(1,(theta - thetar(I))./(thetas(I)-thetar(I)));
            k = KsvG(I) .* S.^LvG(I) .* (1 - (1-S.^(1./mvG(I))).^mvG(I)).^2; 
        end
    end

    function f = aCrop(h)
        % root water uptake reduction factor accordign to Fedds (1978)
        f=ones(size(h));
        f(h> h4 & h<h3) = (h4 - h(h>h4 & h<h3))/(h4-h3);
        f(h<=h4) = 0.0;
    end
    function hacc = hAcc(h)
        % dh/dt for topcel, h=ht(1:2) to update system coefficients of
        % top cell
        [cz,k05,k] = Cond_z(h);
        ccau     = Ccauchy(h(1));
        hacc = (P(it) + cz * (h(2)-h(1)) + ccau(1) * (hCauchy(1)-h(1)) -k05)/   ...
            (dz(1) * Ss(h(1)) + Cpond(h(1)));
    end
    function [cz,K05,k_] = Cond_z(h)
        % Compute or update vertical conductances
        % For top cell use ht(1:2) as input
        N     = numel(h);
        theta = h2theta(h);
        k_     = K(theta);
        
        % implement upwind scheme
        phi   = h+zm(1:N);
        if upWind
            upward       = find( diff(phi) > 0 ); % use head phi = ht+zm 
            downward     = find( diff(phi) < 0 );
            k_(upward)    = k_(upward  +1);
            k_(downward+1)= k_(downward); 
        end
        
        % average k between cell centers
        K05 = (dz(1:N-1)+dz(2:N))./(dz(1:N-1)./k_(1:N-1)+dz(2:N)./k_(2:N)); % harmonic mean

        % conductances
        RZ = 0.5 * dz(1:N) ./ k_;            
        cz = 1 ./ (RZ(1:N-1) + RZ(2:N)); 
    end
    function Cc = Ccauchy(ht)
        %% Set resistance boundary condition for all cells below hDitch

        N = numel(ht);
        In  = find(ht> hCauchy(1:N)); % exfiltrating cells
        Out = find(ht<=hCauchy(1:N)); % infiltrating cells
        Cc= zeros(size(ht));

        % Infiltration
        Cc(In)  = dz(In)./ (Lin.^2./(12.*k(In)*kCover) + Lin .* wIn .* dCover./omega); % active cells with fixed heads
        J       = ht<=hCauchy(1:N) & zm(1:N)>hDitch;  % Infiltration above ditch is impossible so:
        Cc(J)   = 0;

        % Exfiltration
        Cc(Out)  = dz(Out)./(Lex.^2./(12*k(Out)*kCover)); % active cells with fixed heads
        J        = ht>hCauchy(1:N) & zm(1:N)<hDitch; % Exfiltration below ditch level, add exit resistance
        Cc( J)   = dz( J)./(Lex.^2./(12*k( J)*kCover) +  Lex .* wEx .* dCover./omega);

        % Surface runoff and drains
        if ht(1)>0
            Cc( 1)      = Cc(1) + 1./gammaRunoff;
        end
        if iDrain<=N && ht(iDrain)>0
            Cc( iDrain) = Cc(iDrain) + 1./gammaDrain;
        end
    end

    function h = RungeKutta(h0,dt)
        % recursieve Runge-Kutta integratie, 3rd order, O(h^4)
        % Call 3rd order RK function to estimate the next point
        % Verifies that the direct esitmation equals the one done in two
        % stpes of half the time step length each. If not it calls itself
        % with halving the step. This recursive function automatically
        % subdivides the distance to the extent necesassary to attain the
        % required accuracy
        % TO 150928
        dh_ = 0.01;

        h1_ = RK3(h0,dt);           %plot([t,t+dt],[h0,h1],'ro-');
        h2_ = RK3(h0,dt/2);         %plot([t,t+dt/2],[h0,h2],'bo-');
        h3_ = RK3(h2_,dt/2);    %plot([t+dt/2,t+dt],[h2,h3],'go-');
        Delta = (h3_(1)-h1_(1));
        if abs(Delta)>dh_
            h1_ = RungeKutta(h0,dt/2);
            h2_ = RungeKutta(h1_,dt/2);
            h  = h2_;
        else
            h=h3_;
        end
     end 
    function h = RK3(h0, dt)
        %% 3rd order Runge-Kutta integration (Abramowitz and Stegun (1964) 25.5.8)
        k1 = dt * hAcc(h0);
        k2 = dt * hAcc(h0 + k1/2);
        k3 = dt * hAcc(h0 - k1 + 2*k2);
        h  = [h0(1) + k1/6 + 2*k2/3 + k3/6; h0(2)];
    end
end



classdef solutionObj
    %SOLUTION class definition for objects of anlytical solutions for dynamics
    % of flow in a 1D top layer of relatively low transmissivity resting on
    % a regional aquifer with large transmissivity, so that the boundary
    % condition between the two layers is either the head in the second
    % regional aquifer of the (upward positive and uniform) seepage from
    % the second aquifer into the top layer.
    %
    % Class def for anlytical solution objects to compute dynamics of grw
    % fluctusations in a 1 and 2 layer aquifer
    % Inlcludes one numeric solution for the dynamics of the average head in a one- or twolayer cross section
    % TO 110714 110917 (now works for single cross section)
    properties
        name      % name of the solution
        id        % id of the cross section
        b         % [ m ] half width of parcel
        
        dd        % [ m ] ditch depth
        wd        % [ m ] ditch width
        % omega = dd + wd/2 for half dictch.
        
        GR_ELEV   % [NAP] elevation of ground surface

        t         % [ d ]time
        N         % [m/d] net recharge
        q         % [m/d] seepage (upward positive), vector of equal length of N

        D         % [ m; m ] effective thickness of the two aquifers
        Kx        % [m/d ; m/d] horizontal conductivity in the two aquifers
        Ky        % [m/d ; m/d] vertical   conductivity in the two aquifers
        kD        % [m2/d; m2/d] transmissivity of the two aquifers
        c         % [ d ]  resistance between the two aquifers
        
        S         % [ -; - ]  storage coefficient of the two aquifers
        cex       % [ t ]  exfiltration resistance of ditch bottom and side
        cin       % [ t ]  infiltration resistance of ditch bottom and side
        h_mean    % [NAP]  average ditch level
        h_summer  % [NAP]  summer  ditch level
        h_winter  % [NAP]  winter  ditch level
        Qsw       % surface water runoff
        ht        % transient head, average within cross section
        hx        % steady head on input of last time step, varies along cross section
        he        % mean final head based on last time step
        xGr       % [ m ] coordinates along x-axis, only necessary for hx
        Z         % [NAP; NAP] top and bottom of the 2 aquifes (2 values per section)
        HY        % Hydroloigcal Year. It is a struct with length equal to the
                  % number of hydrological years in the data, containing for yeach
                  % hydrological year the times and values of the lowest, highest
                  % and spring head values (HY(iyr).ghg, ghv and glg)
                  % The GHG, GVG and GLG over the last 8 years of the data
                  % are automatically competed from these values.
                  % Note that the spring GVG values are computed as average
                  % over the values at 14 March, 28 March and 14 april of
                  % the same year, not hydrological year.
        GHG       % [ m ]  average highest groundwater level (last 8 hydrological years)
        GVG       % [ m ]  average spring  groundwater level (last 8 hydrological years)
        GLG       % [ m ]  average lowest  groundwater level (last 8 hydrological years)

    end
    properties (Constant)
        solutions={
            % solutions defines the names of the solutions implemented,
            % both their short name and their long name. To tell the
            % solution object is a so and so solution use its short name.
            'Conv'   'Confined exact by convolution';
            'Cos'    'Confined cos approx';
            'Conf'   'Confined seepage added to recharge';
            'Semi'   'Semiconf on aquitard on constant head specified average seepage';
            'GGOR'   'Semiconf on aquitard on constant head specified average seepage, GGOR Ouboter';
            'GGOR2'  'Semiconf on aquitard with given seepage and ditches potentially penetrating';
            'Semiw'  'Semiconf on aquitard with specified seepage with ditch entry resistance;'                    
            'Numeric'      '2 Layer numerical fdm with entry resistances';
            'Semi2w_ode45' '2 layer solution using numerical ode45 solver';
            'Semi2w_rk3'   '2 layer solution using numericla 3rd order Runge-Kutta method';
            'Semi2w_xi'    '2 layer solution direct integration using intermediate variable xi';
            'Semi2w'       '2 layer solution direct integration';
            };
    end
    properties (Dependent=true)
        xm
        Nx
        Dx  
        lambda    % [ m ]  sqrt(kD(:,1).*c)
        gamma     % [ - ]  b/lambda
        Lambda    % [ - ]  tanh(gamma)/gamma
        Gamma     % [ - ]  Gamma/(1-Gamma)
        win       % [ t ]  entry resistance for infiltrating ditch (only implemented for solution semiw (experimental)
        wex       % [ t ]  exfiltration resistance of ditch (check definition), exfiltration=infiltration=default
    end 

    methods

        function o=solutionObj(P)
            % Constructor
            % o is the object instance, same as self in Python
            % P is a struct array of properties one element per case.
            %
            % P has the following field members:
            %  P.L   =    distance between parallel ditches (width of cross section)'
            %  P.AHN =    ground surface elevation
            %  P.q   =    time average given upwards seepage from regional aquifer
            %  P.D1effective  =    effective thickness of first layer (cover layer)
            %  P.D2  =    thickness of regional aquifer, (second layer)
            %  P.hk1 =    horizontal conductivity of first layer
            %  P.hk2 =    horiozntal conductivity of second layer
            %  P.C   =    hydraulic resistance between cover layer and regional aquifer
            %  P.dd  =    ditch depth relative to ground surface !!
            %  P.wd  =    ditch width
            %  P.MU  =    specific yield of cover layer
            %  P.S2  =    elastic storage coefficient of regional aquifer
            %  P.cex =    exfiltration resistance of ditch bottom
            %  P.cin =    infiltration resistance of ditch bottom
            %  P.GP  =    mean surface water head elevation (stage)
            %  P.ZP  =    mean summer surface water head elevation (1 Apr - 30 Sept)
            %  P.WP  =    mean winter surface water head elevation (1 Oct - 31 March) 

            if nargin==0
                return;        % return empty solution obj
            elseif numel(P)>1, % simultaneously generate a list of solutions
                for i=length(P):-1:1
                    o(i)=solutionObj(P(i));
                end
                return;
            else
                o.name=    '';
                o.id  =    P.FID3; % id of cross section, pick from P
                o.b   =    P.L / 2; % half width of parcel
                o.GR_ELEV =P.AHN; % ground surface elevation
                o.t   =    [];    % input time series
                o.N   =    [];    % input precipiation surplus series
                o.q   =    P.q;   % seepage upward from regional aquifer
                o.D   =    [P.D1effective  P.D2]; % Layer thickness
                o.Kx  =    [P.hk1 P.hk2];  %layer hor. conductivities
                o.kD  =    o.Kx .* o.D;      % transmissivity of cover layer
                o.c   =    P.C;            % resistance between cover layer and regional aqufer
                o.Ky  =    [o.D(1)./o.c, o.Kx(2)] / 2; % Vertical conductivity of cover layer
                o.dd  =    P.dd;  % ditch depth with respect to ground surface
                o.wd  =    P.wd;  % ditch width, considered equal to Omega

                o.S   =    [P.MU P.S2];  % storage coefficients
                o.cex =    P.cex;        % exfiltration resistance of ditch bottom
                o.cin =    P.cin;        % infiltration resistance of ditch bottom
                            
                o.h_mean  = P.GP; % average surface water elevation
                o.h_summer= P.ZP; % average surface water elevation summer
                o.h_winter= P.WP; % average surface water elevation winter
                o.Qsw = [];   % surface water runoff
                o.ht  = [];   % transient mean head
                o.hx  = [];   % steady head on input of last time step
                o.he  = [];   % mean final head based on last time step
                o.xGr = [];   % coordinates of x-axis
                o.Z   = [];   % layer elevation
                o.HY  = [];   % horizontald conductivty of phreatic aquifer (MODFLOW, old module)

                % Elevation of tops and bottoms of layers
                z=repmat(o.h_mean,[1,1,3]);
                z(:,:,2) = z(:,:,1)-o.D(:,1);
                z(:,:,3) = z(:,:,2)-o.D(:,2);
                o.Z=z;
            end
        end
        function o=solve(o,name,tne,xGr)
            % Solves the case and stores the results in the object.
            % USAGE obj = obj.solve(name,tne,xGr)
            %  obj is the solution obj
            %  name is the solution short name as specified in "solutions"
            %  tne is a time series with time, precipitation and
            %  evaporation columns, assuming recharge equals precip - evap.
            % xGr is the x-coordinates for which the solution is desired.
            % The results can be plotted on thse coordinate locations.
            
            if numel(o)>1
                for io=numel(o):-1:1
                    o(io) = o(io).solve(name,tne,xGr);
                end                
                return
            end
            
            o=o.legal(name); % assert that name is one of the implemented names/solutions
            o.name = name;
            
            fprintf('Solving cross section with FID = %d, method is %s\n',o.id,o.name);
            
            % Default grid in case grid was not specified.
            if ~exist('xGr','var')
                dx=1.0;         % hard wired
                xGr=0:dx:o.b;
            end
            o.xGr = unique([0, xGr, o.b]);
            
            Nz    = size(o.Z,3)-1;
                        
            % To get a head value at the end of every day, including the
            % first day, use dummy for Dt(1). Length of vector Dt then
            % equals length of vector t.
            Dt    = [NaN; diff(tne(:,1))];  Dt(1)=Dt(2);
            o.N   = tne(:,2)-tne(:,3);  % recharge
            Nt    = length(Dt);
            o.Qsw = zeros(1,Nt);        % surface runoff

            % summer and winter ditch level
            o.t     = tne(:, 1);          % times (actually datenums)
            [~,MON] = datevec(tne(:, 1)); % get months
            
            o.ht = zeros(Nz, length(tne(:,1))); % Initialize time series
            
            switch lower(o.name)
                case 'conv'
                    % Single aquifer confined by convolution derived from
                    % Kraaijnhof van der Leur, series solution. Exact.
                    fprintf('Warning, %s cannot compute surface runoff!\n',o.name);

                    T    =4*o.b^2*o.S(1)/(pi^2*o.kD(1));

                    n=14; tau  =0:Dt(1):n*T; % derived n=7, taken 2x to be sure

                    hLR = ones(size(tne(:,1))) * o.h_summer; hLR(tne(:,end)==0)= o.h_winter;
                    
                    % here we do time at once, but loop over the cross
                    % sections
                    o.ht=hLR(1)+...
                        filter(BR_recharge(T,o.S(1),tau),1,o.N+o.q)'+...
                        filter(BR_ditch(T,tau),1,hLR(:)-hLR(1))';
                    
                    % steady state
                    F=hLR(end)+3/2*(o.N(end)+o.q(end)).*(o.b.^2./(3*o.kD(:,1)));
                    for ix=1:length(o.xm)
                        o.hx(:,ix)=F.*(1-(o.xm(ix)./o.b).^2);
                    end
                    o.he =((o.N(end)+o.q(end)).*o.b^2)./(3*o.kD(:,1)); 

                    
                case 'cos'  % Single aquifer cosine approx, i.e. first term of
                            % Kraaijenhof van der Leur. Seepage directly added to recharge
                    
                    T    =4/(pi^2)*o.b.^2.*o.S(1)./o.kD(1);

                    for it=1:length(Dt);
                        itq = max(numel(o.q),it);
                        hLR=o.h_winter+(MON(it)>=4 && MON(it)<=9).*(o.h_summer-o.h_winter);
                        e=exp(-Dt(it)./T);
                        if it==1
                            o.ht(:,it)=hLR + (o.N(it)+o.q(itq)).*T./o.S(:,1).*(1-e);
                        else
                            o.ht(:,it)=hLR +(o.ht(:,it-1)-hLR).*e + (o.N(it)+o.q(itq)).*T./o.S(:,1).*(1-e);
                        end
                        % tackle surface water runoff
                        if o.ht(1,it)>o.GR_ELEV
                            o.Qsw(it)=o.b.*o.S.*(o.ht(1,it)-o.GR_ELEV);
                            o.ht(1,it)=o.GR_ELEV;
                        end
                    end
                    
                case 'conf'  % Single confined aquifer seepage directly added to recharge.
                    
                    T   = o.b.^2.*o.S(1)./(3*o.kD(1));  % characteristic time of this solution
                    
                    for it=1:length(Dt);
                        itq = max(numel(o.q),it);
                        hLR=o.h_winter+(MON(it)>=4 && MON(it)<=9).*(o.h_summer-o.h_winter);
                        e=exp(-Dt(it)./T);
                        if it==1
                            o.ht(:,it)=hLR+(o.N(it)+o.q(itq)).*T./o.S(:,1).*(1-e);
                        else
                            o.ht(:,it)=hLR +(o.ht(:,it-1)-hLR).*e + (o.N(it)+o.q(itq)).*T./o.S(:,1).*(1-e);
                        end
                        % tackle surface water runoff
                        if o.ht(1,it)>o.GR_ELEV
                            o.Qsw(it)=o.b.*o.S(1).*(o.ht(1,it)-o.GR_ELEV);
                            o.ht(1,it)=o.GR_ELEV;
                        end
                    end
                    
                case 'semi'  % Aquif on top of aquitard below which constant head

                    T      =o.b.^2.*o.S(1)./(3*o.kD(1)).*o.Gamma;
                    
                    for it=1:length(Dt);
                        itq = max(numel(o.q),it);
                        hLR=o.h_winter+(MON(it)>=4 && MON(it)<=9).*(o.h_summer-o.h_winter);
                        e=exp(-Dt(it)./T);
                        if it==1
                            o.ht(:,it) =hLR +(o.N(it)+o.q(itq)).*T./o.S(:,1).*(1-e);
                        else
                            o.ht(:,it)=hLR +(o.ht(:,it-1)-hLR).*e + (o.N(it)+o.q(itq)).*T./o.S(:,1).*(1-e);
                        end
                        % tackle surface water runoff
                        if o.ht(1,it)>o.GR_ELEV
                            o.Qsw(it)=o.b.*o.S(1).*(o.ht(1,it)-o.GR_ELEV);
                            o.ht(1,it)=o.GR_ELEV;
                        end
                    end
                    
                    % steady state
                    R=(o.N(end)+o.q(end)).*o.c;
                    P=cosh(o.gamma);
                    for ix=1:length(o.xm)
                        o.hx(ix)   =R.*(1-cosh(o.xm(ix)./o.lambda)./P);
                    end
                    o.he   =hLR+R.*o.Gamma;
                    
                case 'ggor' % single aquifer on top of aquitard below which constant head
                            % As applied in the GGOR-tool bij Waternet (till 2011)
                            % TO 101113 101114 110106

                    EPS=0.9; EPS_m1=1/EPS; % implicitness
                    F = o.Lambda./(1-o.Lambda);
                    
                    for it=1:length(Dt)
                        itq = max(numel(o.q),it);
                        hLR=o.h_winter+(MON(it)>=4 && MON(it)<=9).*(o.h_summer-o.h_winter);
                        Fmuc=o.c.*o.S(:,1)/Dt(it);
                        if it==1,
                            o.ht(:,it)=EPS_m1*(hLR.*Fmuc+hLR.*F+...
                                o.c.*(o.N(it)+o.q(itq)))./(F+Fmuc)+(1-EPS_m1)*hLR;
                        else
                            o.ht(:,it)=EPS_m1*(o.ht(:,it-1).*Fmuc+hLR.*F+...
                                o.c.*(o.N(it)+o.q(itq)) )./(F+Fmuc)+...
                                    (1-EPS_m1)*o.ht(:,it-1);
                        end
                        % tackle surface water runoff
                        if o.ht(1,it)>o.GR_ELEV
                            o.Qsw(it)=o.b.*o.S(1).*(o.ht(1,it)-o.GR_ELEV);
                            o.ht(1,it)=o.GR_ELEV;
                        end                       
                    end
                    
                    % steady state
                    R=(o.N(end)+o.q(end)).*o.c;
                    P=cosh(o.gamma);
                    for ix=1:length(o.xm)
                        o.hx(ix)   =R.*(1-cosh(o.xm(ix)./o.lambda)./P);
                    end
                    o.he   =hLR+R.*o.Gamma;

           
                case 'semiw';                    
                    % Single aquifer on top of aquitad below which constant head
                    % seepage from below is specified
                    % The ditches have entry resistance (same as the multi-layer solution)
                    % however, in case the ditch infiltrates a different
                    % resistance is used, yielding B1 and T1 instead of B
                    % and T

                    B = (1 - ...
                        (o.b./o.lambda*(o.Kx(1).*o.wex(1)./o.lambda+...
                        1./tanh(o.b./o.lambda)))...
                        .^(-1)).^(-1);
                    T = o.c .* o.S(1) * 1 ./ (B -1);
                    
                    % in case of infiltration from ditch use different
                    % ditch bottom resistance win
                    B1= (1 - ...
                        (o.b./o.lambda*(o.Kx(1).*o.win(1)./o.lambda+...
                        1./tanh(o.b./o.lambda)))...
                        .^(-1)).^(-1);
                    T1= o.c .* o.S(1) * 1 ./ (B1-1);

                    for it=1:length(Dt) 
                        itq = max(numel(o.q),it);
                        hLR=o.h_winter+(MON(it)>=4 && MON(it)<=9).*(o.h_summer-o.h_winter);
                        
                        % assume flow from aquifer to ditch
                        e=exp(-Dt(it)./T);
                        if it==1,
                            o.ht(:, 1)=hLR+(o.N(it)+o.q(itq)).*T./o.S(:,1).*(1-e);
                        else
                            o.ht(:,it)=hLR+(o.ht(:,it-1)-hLR).*e + (o.N(it)+o.q(itq)).*T./o.S(:,1).*(1-e); 
                        end
                        
                        if o.ht(1,it)<hLR % in case flow from ditch into aquifer
                            e=exp(-Dt(it)./T1);
                            if it==1,
                                o.ht(:, 1)=hLR+(o.N(it)+o.q(itq)).*T1./o.S(:,1).*(1-e);
                            else
                                o.ht(:,it)=hLR+(o.ht(:,it-1)-hLR).*e + (o.N(it)+o.q(itq)).*T1./o.S(:,1).*(1-e); 
                            end
                        end

                        
                        % tackle surface water runoff
                        if o.ht(1,it)>o.GR_ELEV;
                            o.Qsw(it)=o.b*o.S(1).*(o.ht(1,it)-o.GR_ELEV);
                            o.ht(1,it)=o.GR_ELEV;
                        end
                    end
                    
                   % steady state
                    R=(o.N(end)+o.q(end)).*o.c;
                    P=o.Kx(:,1).*o.win(:,1)./o.lambda.*sinh(o.gamma)+cosh(o.gamma);
                    for ix=1:length(o.xm)
                        o.hx(:,ix)   =R.*(1-cosh(o.xm(ix)./o.lambda)./P);
                    end
                    o.he   =hLR+R.*T./(B-1);
 
               case 'GGOR2';
                    % As semiw, but with ditches potentially penetrating
                    % into the second aquifer. Flow from second aquifer to
                    % the ditches is parallel to flow through aqutard.
                    % Ditches have in- and exfiltration resistance.
                    % Single aquifer solution with total seepage specified.
                    % The total seepage is the sum of the seepage through
                    % the aquitard and that directly from the second
                    % aquifer to penetrating ditches. The entry resistance
                    % of the ditches in the second aquifer is such that it
                    % also works when ditches are not penetrating, which
                    % yields an entry resistance of the ditches in the
                    % second aquifer equal to infinitiy.

                    F = o.wex(1)/o.b + o.b/o.lambda * coth(o.b/o.lambda)-1;
                    F1= o.win(1)/o.b + o.b/o.lambda * coth(o.b/o.lambda)-1;
                    
                    B = (1 - ...
                        (o.b./o.lambda*(o.Kx(1).*o.wex(1)./o.lambda+...
                        1./tanh(o.b./o.lambda)))...
                        .^(-1)).^(-1);
                    T = o.c .* o.S(1) * 1 ./ (B -1);
                    
                    % in case of infiltration from ditch use different
                    % ditch bottom resistance win
                    B1= (1 - ...
                        (o.b./o.lambda*(o.Kx(1).*o.win(1)./o.lambda+...
                        1./tanh(o.b./o.lambda)))...
                        .^(-1)).^(-1);
                    T1= o.c .* o.S(1) * 1 ./ (B1-1);
                    
                    
                    for it=1:length(Dt) 
                        itq = max(numel(o.q),it);
                        hLR=o.h_winter+(MON(it)>=4 && MON(it)<=9).*(o.h_summer-o.h_winter);
                        
                        fac = o.c/o.b*o.wex(2)/o.D(2);
                        
                        % assume flow from aquifer to ditch
                        e=exp(-Dt(it)./T);
                        if it==1,
                            o.ht(:, 1)=hLR+(o.N(it)+o.q(itq)/(1+fac)).*T./o.S(:,1).*(1-e);
                        else
                            o.ht(:,it)=hLR+(o.ht(:,it-1)-hLR).*e + (o.N(it)+o.q(itq)/(1+fac)).*T./o.S(:,1).*(1-e); 
                        end
                        
                        if o.ht(1,it)<hLR % in case flow from ditch into aquifer
                            fac= o.c/o.b * o.win(2)/o.D(2);
                            
                            e=exp(-Dt(it)./T1);
                            if it==1,
                                o.ht(:, 1)=hLR+(o.N(it)+o.q(itq)/(1+fac)).*T1./o.S(:,1).*(1-e);
                            else
                                o.ht(:,it)=hLR+(o.ht(:,it-1)-hLR).*e + (o.N(it)+o.q(itq)/(1+fac)).*T1./o.S(:,1).*(1-e); 
                            end
                        end

                        
                        % tackle surface water runoff
                        if o.ht(1,it)>o.GR_ELEV;
                            o.Qsw(it)=o.b*o.S(1).*(o.ht(1,it)-o.GR_ELEV);
                            o.ht(1,it)=o.GR_ELEV;
                        end
                    end
                    
                   % steady state
                    R=(o.N(end)+o.q(end)).*o.c;
                    P=o.Kx(:,1).*o.win(:,1)./o.lambda.*sinh(o.gamma)+cosh(o.gamma);
                    for ix=1:length(o.xm)
                        o.hx(:,ix)   =R.*(1-cosh(o.xm(ix)./o.lambda)./P);
                    end
                    o.he   =hLR+R.*T./(B-1);
 
                case 'numeric' % Numerical 2Layer solution with entry resistances

                    % problem tackled on a per cross section basis to allow
                    % using the two-D matlab finite difference models

                    o.ht = zeros(Nz,Nt);                                        
%                     o.hx = zeros(Nz,length(o.xm));
%                     o.he = zeros(Nz,1);
                    
                    FQ = zeros(Nz,Nx);
                    FH = NaN(Nz,Nx);
 
                    zGr = squeeze(o.Z(1,:));

                    kx = o.Kx(:)*ones(size(Dx));
                    C  = o.c    *ones(size(Dx));
                    ky = Inf;
                    
                    wequiv = o.win(:)+0.5*Dx(1)./kx(:,1);

                    kx(:,1)=0.5*Dx(1)./wequiv(:);

                    Ss  = (o.S(:)./o.D(:))*ones(size(Dx));

                    %% Initial and boundary conditions

                    for it = 1:length(Dt)
                        itq = max(numel(o.q),it);
                        hLR = o.h_winter+(MON(it)>=4 & MON(it)<=9)*(o.h_summer-o.h_winter);
                        FQ(1,:) = o.N(it).*Dx;
                        FQ(2,:) = o.q(itq)    .*Dx;

                        if it==1
                            IH = ones(Nz,Nx) *hLR;
                            FH(:,1) = IH(:,1);
                            FI = fdm2ct(o.xGr,zGr,o.t(1:2),kx,C,ky,Ss,IH,FH,FQ); % compute end of first day separately
                        else
                            FH(:,1) = hLR;
                            FI = fdm2ct(o.xGr,zGr,o.t([it-1,it]),kx,C,ky,Ss,FI(:,:,end),FH,FQ);
                        end

                        % tackle surface runoff
                        fi = FI(1,:,end);
                        I  = find(fi>o.GR_ELEV);
                        if ~isempty(I)
                            o.Qsw(it) = sum((fi(I)-o.GR_ELEV).*((Ss(1,I)*Obj.D(1)).*Dx(I)));
                            fi(I) = o.GR_ELEV;                            
                            FI(1,:,end) = fi;
                        end                           
                        o.ht(1,it) = sum(FI(1,:,end).*Dx,2)/o.b;
                        o.ht(2,it) = sum(FI(2,:,end).*Dx,2)/o.b;
                        
                        % Steady state
                        hx_ = squeeze(FI(:,:,end));
                        o.he(1) = sum(hx_(1,:).*Dx)/o.b;
                        o.he(2) = sum(hx_(2,:).*Dx)/o.b;

                    end
                    fprintf('\n');
                    
                    % simulation results for cross section
                    o.hx = fdm2c(o.xGr,zGr,kx,C,ky,FH,FQ);
                    o.he = sum(o.hx.*[Dx; Dx],2)./o.b;

                otherwise
                    % Analytic solution for single aquifer on top of seminconfined aquifer separated
                    % by leaking semi-confined bed with given supply in second aquifer.
                    %
                    % TO 101113 101114 101211
                    
                    NLay=2; I = eye(NLay);
                    o.ht=zeros(NLay,length(Dt));
                    o.hx=zeros(NLay,length(o.xm));
                    o.he=zeros(NLay,1);

                        
                        hLR = o.h_winter+(MON>= 4 & MON<=9)*(o.h_summer-o.h_winter);
                        kd  = o.kD(:);
                        cc   = [1e6 o.c Inf]';
                        s    = o.S(:);
                        
                        T    = diag(kd); %T_m1 = T^(-1);
                        s    = diag(s);  S_m1 = s^(-1);
                        H_m1 = diag(o.Kx(:)'.*o.win(:)');

                        A = -diag( 1./(kd(2:NLay  ).*cc(2:NLay)),-1)+...
                            +diag( 1./(kd(1:NLay  ).*cc(1:NLay))+1./(kd(1:NLay).*cc(2:NLay+1)), 0)+...
                            -diag( 1./(kd(1:NLay-1).*cc(2:NLay)),1);

                        % A_m1  = A^(-1);
                        sqrtA = sqrtm(A); sqrtA_m1=sqrtA^(-1); % sqrtA_m1= sqrtA^(-1);

                        sinhm = (expm(o.b*sqrtA)-expm(-o.b*sqrtA))/2;
                        coshm = (expm(o.b*sqrtA)+expm(-o.b*sqrtA))/2;

                        F    = (H_m1*sqrtA*sinhm+coshm);
                        F_m1 = F^(-1);
                        TAB    = T*A*(I-sqrtA_m1/o.b * sinhm * F_m1)^(-1);
                        G = S_m1*TAB;
                        [V,E] = eig(G);
                        V_m1 = V^(-1);
                        E_m1 = E^(-1);
                        
                        h = zeros(2,Nt);
                        
                        switch lower(o.name)
                            case 'semi2w_ode45'
                                fprintf('Warning, Semi2w_ode5 cannot compute surface runoff!\n');

                                [~,h] = ode45('solveNLay',tne(:,1),[hLR(1); hLR(1)],odeset,S_m1,TAB,tne,o.q,hLR);
                                o.ht(:) = h';
                                
                            case 'semi2w_rk3'  %% Diect intration by third order Runge Kutta method with step control  
                                
                                for it = 1:length(Dt)
                                    itq = max(numel(o.q),it);
                                    fprintf('.');
                                    hLR = o.h_winter+(MON(it)>=4 && MON<=9)*(o.h_summer-o.h_winter);
                                    if it==1
                                        H0 =[hLR;hLR];
                                    else
                                        H0 = h(:,it-1);
                                    end
                                    dt = Dt(it); tau = 0;
                                    while tau<Dt(it)
                                        f1 = S_m1*([o.N(it);o.q(itq)]-TAB*(H0-hLR)); H1=H0+f1*dt/2;
                                        f2 = S_m1*([o.N(it);o.q(itq)]-TAB*(H1-hLR)); H2=H0+dt*(2*f2-f1);
                                        f3 = S_m1*([o.N(it);o.q(itq)]-TAB*(H2-hLR));
                                        H3 = H0+dt*(f1+4*f2+f3)/6;
                                        d  = dt*(f1-2*f2+f3)/6; 
                                        if  max(abs(d))<1e-4
                                            tau = tau+dt;
                                            H0 = H3;
                                        else
                                            dt = dt/2;
                                        end
                                    end
                                    h(:,it) = H3;
                                    % tackle surface water runoff
                                    if h(1,it) > o.GR_ELEV
                                        o.Qsw(it) = o.b*o.S(1)*(h(1,it)-o.GR_ELEV);
                                        h(1,it) = o.GR_ELEV;
                                    end

                                    if rem(it,100) == 0, fprintf('%d\n',it); end
                                end
                                
                            case 'semi2w_xi' % Making the equations independent by means of eigen vectors and eigen values
 
                                xi     =zeros(2,Nt);
                                for it  = 1:length(Dt)
                                    itq = max(numel(o.q),it);
                                    hLR = o.h_winter+(MON(it)>=4 && MON(it)<=9)*(o.h_summer-o.h_winter);
                                    e = expm(-E*Dt(it));
                                    theta = V_m1*S_m1*[o.N(it);o.q(itq)];
                                    if it == 1
                                        xi(:,it) = (I-e)*E_m1*theta;
                                    else
                                        xi(:,it) = e*xi(:,it-1)+(I-e)*E_m1*theta;
                                    end
                                    h(:,it) = hLR+V*xi(:,it);
                                    % tackle surface water runoff
                                    if h(1,it) > o.GR_ELEV
                                        o.Qsw(it) = o.b*o.S(1)*(h(1,it)-o.GR_ELEV);
                                        h(1,it)   = o.GR_ELEV;
                                        xi(:,it)  = V_m1*(h(:,it)-hLR);
                                    end
                                end

                            case  'semi2w'; % omit the step via xi
                                R = E_m1*V_m1*S_m1;
                                
                                for it = 1:length(Dt)
                                    itq = max(numel(o.q),it);
                                    hLR = o.h_winter+(MON(it)>=4 && MON(it)<=9)*(o.h_summer-o.h_winter);
                                    e = expm(-E*Dt(it));
                                    if it==1
                                       h(:,it) = hLR+V*(I-e)*R*[o.N(it);o.q(itq)];
                                    else
                                       h(:,it) = hLR+V*e*V_m1*(h(:,it-1)-hLR)+V*(I-e)*R*[o.N(it);o.q(itq)];
                                    end
                                    % tackle surface water runoff
                                    if h(1,it)>o.GR_ELEV
                                        o.Qsw(it)=o.b*o.S(1)*(h(1,it)-o.GR_ELEV);
                                        h(1,it)=o.GR_ELEV;
                                    end
                                end
                            otherwise
                                error('un recognized solution method requested %s, use one of ...\n%s',...
                                    o.name,sprintf(' ''%s''',o.solutions{:,1}));
                        end
                        o.ht(1,:)=h(1,:);
                        o.ht(2,:)=h(2,:);
                        
                        % steady state
%                          hx_=zeros(NLay,length(o.xm));
%                          N_ = mean(o.N);
%                          q_ = mean(o.q);
%                          for ix=1:length(o.xm);
%                              hx_(:,ix)=hLR+(I-funm(o.xm(ix)*sqrtA,@cosh)*F_m1)*A_m1*T_m1*[N_;q_];
%                          end
%                          o.hx=permute(hx_,[3,2,1]);
%                          o.he=sum(hx_.*(ones(Nz,1)*Dx),2)./o.b;
            end
             
            [o.GLG,o.GVG,o.GHG,o.HY] = getGXG(o.ht(1,:),o.t);
            
        end
        
%         function o=getGXG(o)
%             %[GLG,GVG,GHG]=getGXG(h,t,plotmode)
%             % Compute Average lowest head (GLG), spring head (GVG) and highest head
%             % (GHG) from the given heads over the given time span.
%             % the heads are presumed ordere in the direction of the time vector
%             % therefore, a large number of time series can be handled at once
%             % plot=0 or omitted altogether, don't plot
%             % plot=1 plot heads over time with points that were used to compute CxG
%             % plot=2 plot GxG for all sections with section number on x-axis
%             % TO 101023
% 
%             if numel(o)>1
%                 for io=numel(o):-1:1
%                     o(io)=getGXG(o(io));
%                 end
%                 return
%             end
%             
%             DV=datevec(o.t(:));
% 
%             % get hydrological years in time series
%             if o.t(  1)<=datenum(DV(  1,1),4, 1),  % if first date < April 1 in year 1
%                 yr1 = DV(  1,1);
%             else
%                 yr1 = DV(  1,1)+1;
%             end
%             if o.t(end) >= datenum(DV(end,1),3,31), % if last data > 31 March last year
%                 yr2 = DV(end,1);
%             else
%                 yr2 = DV(end,1)-1;
%             end
% 
%             yrs = yr1:yr2;  % hydrological years in time series
%             
%             for ihy = length(o.HY):-1:1    % for all hydrological years in data
% 
%                 % GVG
%                 K=find( ...
%                         DV(:,1)==yrs(ihy) & (...
%                            ( DV(:,2)==3 & ( DV(:,3)==14 | DV(:,3)==28 ) ) |...
%                            ( DV(:,3)==4 &   DV(:,3)==1 ) ...
%                                             ) ...
%                       );
%                 o.HY(ihy).gvg = [o.t(K), o.ht(K)];
% 
%                 % GHG and GLG
%                 o.HY(ihy).t1=datenum(yrs(ihy)  ,4,1);  % its starting datenum
%                 o.HY(ihy).t2=datenum(yrs(ihy)+1,4,1);  % its ending   datenum
%                 o.HY(ihy).J=find( o.t>=o.HY(ihy).t1 &...
%                                    o.t<o.HY(ihy).t2 &...
%                                   (DV(:,3)==14 | DV(:,3)==28));
% 
%                 R=sortrows([o.t(o.HY(ihy).J) o.ht(1,o.HY(ihy).J)'],2);
% 
%                 o.HY(ihy).glg=R(1:3      ,:);
%                 o.HY(ihy).ghg=R(end-2:end,:);
%             end
%         end      
        
        function o = showGXG(o)
        %% Show data and points picked out for GXG computation
            yrs=datevec([o(1).HY.t1]);

            figure; hold on;
            xlabel(sprintf('%d - %d',yrs(1),yrs(end)));
            ylabel('head [m]'); grid on;
            title(sprintf('Solution=%s: head and values picked out for GXG',o(1).name));
            for io=numel(o):-1:1
                plot(o(io).t([1 end]),[o(io).GHG o(io).GHG],'r--');
                plot(o(io).t([1 end]),[o(io).GVG o(io).GVG],'g--');
                plot(o(io).t([1 end]),[o(io).GLG o(io).GLG],'b--');

                J=vertcat(o(io).HY.J);
                plot(o(io).t(J),o(io).ht(1,J), 'bo'); % plot all

                plot(o(io).GHG(:,1),o(io).GHG(:,2),'ro','markerFaceColor','r');
                plot(o(io).GVG(:,1),o(io).GVG(:,2),'go','markerFaceColor','g');
                plot(o(io).GLG(:,1),o(io).GLG(:,2),'bo','markerFaceColor','b');
                plot(o(io).t,o(io).ht(1,:),'k');
            end
            datetick('x',12);

            legend('GHG','GVG','GLG','sim14','GHGyr','GVGyr','GLGyr','sim');
        end

        function o=plott(o,clr,linestyle,t)
            if numel(o)>1
                ttl=o(1).name;
                title(sprintf('Solution %s',ttl));
                xlabel('time'); ylabel('head [m]');
                hold on; grid on;
                for io=1:numel(o)
                    o(io)=plott(o(io),clr,linestyle,t);
                end
                return;
            end
            
            twolayers = size(o.ht,1)>1;
            if nargin<4, t=o.t; end
            if nargin<3, linestyle='-'; end
            if nargin<2, clr='r'; end
            plot(t,o.ht(1,:),[clr(1) linestyle(1)]);
            if twolayers
                plot(t,o.ht(2,:),[clr(end) linestyle(end)])
            end
            %datetick;
        end
        
        function o=plotx(o,clr)
            hold on; grid on;
            title(sprintf('Solution %s, steady XSec',o(1).name));
            xlabel('x [m]'); ylabel('head [m]');

            if numel(o)>1
                for io=1:numel(o)
                    o(io)=plot(o(io),clr);
                end
                return;
            end
            twolayers = length(size(o.hx(:,1)))>1;
            if nargin==1, clr='br'; end
            plot(o.xm,o.hx(1,:),clr(1));
            plot(o.xGr([1 end]),o.he(1,[1 1]),[clr(1) '--']);
            if twolayers
                plot(o.xm,o.hx(2,:),clr(end));
                plot(o.xGr([1 end]),o.he(2,[1 1]),[clr(end) '--']);
            end
        end
        
        function br=BR_recharge(T,mu,tau)
            % BR_recharge block response for recharge
            NRech=1;
            n=max(30,1.5*T/diff(tau(1:2))); % afgeleide criterium bleek te zwak, 30 door experiementeren
            sr=(1-exp(-tau/T));
            for j=2:n
                j2=2*j-1;
                sr=sr+(1/j2)^4*(1-exp(-j2^2*tau/T));
            end
            sr=NRech*T/mu*8/(pi^2)*sr;
            br=diff(sr);
        end
        
        function br=BR_ditch(T,tau)
            % BR_ditch block response for ditch
            A=1;
            n=max(30,1.5*T/diff(tau(1:2))); % afgeleide criterium bleek te zwak, 30 door experimenteren
            sr=exp(-tau/T);
            for j=2:n
                j2=2*j-1;
                sr=sr+(1/j2)^2*exp(-j2^2*tau/T);
            end
            sr=A*(1-8/(pi^2)*sr);
            br=diff(sr);            
        end
        
        function o=legal(o,name)
            if numel(o)>1
                for io=numel(o):-1:1
                    o(io)=legal(o(io),name);
                end
                return;
            end
            if nargin==1
                fprintf('Legal names for the implemented solutions:\n');
                for j=1:numel(o.solutions)
                    fprintf('%15s : %s\n',o.solutions{j,:});
                end
            else % case independent
                I=find(strcmpi(name,o.solutions(:,1)));
                if ~isempty(I)
                    solutionName = o.solutions{I(1),1};
                    o.name = solutionName;
                else
                    error('Solution name <<%s>> not known\n',name);
                end
            end
        end
        
        % Dependent variables/members
        function xm     = get.xm(o), xm = 0.5*(o.xGr(1:end-1)+o.xGr(2:end)); end
        function Dx     = get.Dx(o), Dx = abs(diff(o.xGr)); end
        function Nx     = get.Nx(o), Nx  = numel(o.Dx); end  
        function lambda = get.lambda(o), lambda = sqrt(o.kD(:,1).*o.c); end
        function gamma  = get.gamma(o), gamma =  o.b./o.lambda; end
        function Lambda = get.Lambda(o), Lambda = tanh(o.gamma)./(o.gamma); end % see theory
        function Gamma  = get.Gamma(o), Gamma  = (3./o.gamma.^2).*(1-o.Lambda)./o.Lambda; end
        function wex    = get.wex(o)        % exfiltration resistance of ditches, see theory for background
                 delta  = max(0, o.D(1)-o.dd);
                 lsb1   = min(o.wd, delta);
                 lsb2   = o.wd - lsb1;
                 k1     = sqrt(o.Ky(1)*o.Kx(1));
                 k2     = sqrt(o.Ky(2)*o.Kx(2));
                
                 wex    = [2 * o.D(1) / (pi*k1) * log(pi/2 * o.D(1)/(o.dd+lsb1)) + ...
                          o.D(1) / (lsb1/o.cex + o.dd/o.cex); ...                    
                          2 * o.D(2) / (pi*k2) * log(pi/2 * o.D(2)/(    lsb2)) + ...
                          o.D(1) / (lsb2 * (delta/o.Ky(1) +o.cex)) ...
                       ];
        end
        
        function win = get.win(o) % infiltraiton resistance of ditches, see theory for background
                delta  = max(0, o.D(1)-o.dd);
                lsb1   = min(o.wd, delta);
                lsb2   = o.wd - lsb1;
                k1     = sqrt(o.Ky(1)*o.Kx(1));
                k2     = sqrt(o.Ky(2)*o.Kx(2));
                
                win = [2 * o.D(1) / (pi*k1) * log(pi/2 * o.D(1)/(o.dd+lsb1)) + ...
                       o.D(1) / (lsb1/o.cin + o.dd/o.cin); ...                  
                       2 * o.D(2) / (pi*k2) * log(pi/2 * o.D(2)/(    lsb2)) + ...
                       o.D(1) / lsb2 * (delta/o.Ky(1) +o.cin) ...
                       ];
        end
    end
end


function [H,B]=mFlow(gr,t,IBOUND,kx,ky,kz,Ss,STRTHD,FQ,varargin)
%MFLOW MODFLOW alike 3D block-centred transient finite difference model as a Matlab function for use in class
%
% This model supersedes all fdm models mentioned under See Also below.
%
% Example:
%    [H,B]=mFlow(gr,t,IBOUND,kx,ky,kz,Ss,STRTHD,FQ,varargin)
% 
% INPUT:
%    gr  = gridObj (obtained by calling gr = gridObj(xGr,yGr,zGr)
%          with xGr,yGr,zGr the coordiantes of the grid lines in the 3 dimensions
%          of the model and where zGr may be a full 3D array of valus with
%          Ny and Nx corresponding to cells and Nz corresponding to cell
%          cell interfaces, i.e. size(zGr,3) must be Nz+1);
%    t   = time series (ends of stress periods).
%          if length(t)==0 --> steady state else transient
%          There are no separate options to set steady state or transient,
%          it follow from length(t).
%    IBOUND [ - ] boundary array where IBOUND<0 indicats fixed head, IBOUND>0
%         indicaes head to be computed. IBOUND==0 indicates inactive cells.
%    kx  = [L/T] horizontal conductivity
%    ky  = [L/T] horizontal conductivity in y-direction
%    kz  = [L/T] vertical hydraulic conductivity
%    Ss  = [1/L] specific storage coefficient
%    STRTHD = [L] is a 3D numerical array with the initial heads. These initial 
%         heads correspond with the fixed heads where IBOUND<0.
%    FQ =  [L3/T] is a 3D array with the fixed flows. It allows setting
%         fixed flows for wells and recharge where no further values are
%         necesssary in subssequent time step.
%
%    FUTURE extensions:
%       allow varargin to accept parametername,value pairs to state package
%       information like in MODFLOW, such as for WEL, DRN, GHB, RIV, CHD.
%    ALREADY implemented name value pair to set theta (implicitness)
%       'theta',value -- the implicitness, value between 0.5-1 where 0.67
%
%    Possibilities for input name,value pairs:
%         CHD,values  (MODFLOW-equivalent, to be implmennted)
%         WEL,values  (MODFLOW-equivalent, to be implemented)
%         GHB,values  (MODFLOW-equivalent, to be implemented)
%         DRN,values  (MODFLOW-equivalent, to be implemented)
%         RIV,values  (MODFLOW-equivalent, to be implemented)
%
% OUTPUT
%  H = array of head structs, one element per time step and computed heads
%         stored in H(it).values
%  B = cell by cell flow struct array, one element per time step with
%         cell by cell flows stored in B(it).term{j} where j corresponds to
%         one of the labels in B(it).label. These labels can be one of the
%         labels used in MODFLOW to indicate which flow term is meant:
%            'FLOW RIGHT FACE'
%            'FLOW FRONt FACE'
%            'FLOW LOWER FACE'
%            'CONSsTANt HEAD',
%            'STORAGE'
%            'WELLS'
%
% See also: fmd2 fdm2c fdm2t fdm2NC fdm2tNC fdm2ct fdm2dens fdm3 fdm3t
%
% TO 991017  TO 000530 001026 070414 080301 100307
% TO 080226 (added inactive cells and true fixed heads
% TO 120420  full 3D transient and MOFLOW like, also its outputs

%% Without input, run selftest

if nargin<2
    if nargin==1, i=gr; else i=[]; end
    [H,B] = selfTest(i); % in this case grid is an index 
    return;
end

if ~strcmpi(class(gr),'gridObj')
    error('%s: first argument must be a gridObj',fmilename);
end

%% Modflow compatible labels to indicate cell by cell flow flow type
labels = {'FLOW RIGHT FACE','FLOW FRONT FACE','FLOW LOWER FACE','CONSTANT HEAD','STORAGE','WELLS'};

% column numbers to be used in budget struct B
iRIGHT = strmatchi('FLOW RIGHT',labels);
iFRONT = strmatchi('FLOW FRONT',labels);
iLOWER = strmatchi('FLOW LOWER',labels);
iCONST = strmatchi('CONSTANT'  ,labels);
iSTOR  = strmatchi('STORAGE'   ,labels);
iWELLS = strmatchi('WELLS'     ,labels);

%% implicitess, optional input
[theta,varargin] = getProp(varargin,'theta',0.67);
theta = min(1,max(theta,0.5));

%% get time and see wheather steaddy-state or transient computations are required
if numel(t)>1    % transient flow
    dt=diff(t);
else             % steady-state flow
    dt    = t;
    Ss    = gr.const(0);
    theta = 1;
end

if ~all(dt>0), error('%s: time must be increasing',mfilename); end
Nt = numel(dt);

%% Assert that all these variables are truly 3D arrays
check(gr,'IBOUND',IBOUND,'kx',kx,'ky',ky,'kz',kz,'Ss',Ss,'STRTHD',STRTHD,'FQ',FQ);
        
%% Compute matix coefficients

% Step 1, compute cell resistances
if gr.AXIAL
    RX1 = 1/(2*pi)*log(gr.XM./gr.XGR(1:end-1,1:end-1,1:end-1))./gr.DZ./kx;
    RX2 = 1/(2*pi)*log(gr.XGR(1:end-1,2:end,1:end-1)  ./gr.XM)./gr.DZ./kx;
    RY  = Inf(gr.size);
    RZ  = 0.5*gr.DZ/pi./(gr.XGR(1:end-1,2:end,1:end-1).^2-gr.XGR(1:end-1,1:end-1,1:end-1).^2)./kz;

    Cx=1./(RX1(:,2:end,:)+RX2(:,1:end-1,:));    
    Cy=1./(RY( 1:end-1,:,:)+RY( 2:end,:,:));
    Cz=1./(RZ( :,:,1:end-1)+RZ( :,:,2:end));
else
    RX  = 0.5*(gr.DX./(gr.DY.*gr.DZ))./kx;
    RY  = 0.5*(gr.DY./(gr.DZ.*gr.DX))./ky;
    RZ  = 0.5*(gr.DZ./(gr.DX.*gr.DY))./kz;

    Cx=1./(RX( :,1:end-1,:)+RX( :,2:end,:));    
    Cy=1./(RY( 1:end-1,:,:)+RY( 2:end,:,:));
    Cz=1./(RZ( :,:,1:end-1)+RZ( :,:,2:end));
end

% Setp 2, compute storage matrix coefficient
Cs=gr.Vlay.*Ss;

%% Referencing neighbors of every cell

% step 1 generate cell numbers in the internal order of Matlab arrays
Nodes    = ones(gr.size);    % all 1
Nodes(:) = cumsum(Nodes(:)); % add them in order of natural array order

% get cell nubers of eastern, western, northern, southern, top and bottom neighbors
IE=Nodes(:,2:end,:);   IW=Nodes(:,1:end-1,:);
IS=Nodes(2:end,:,:);   IN=Nodes(1:end-1,:,:);
IT=Nodes(:,:,2:end);   IB=Nodes(:,:,1:end-1);

%% Put the matrix coefficients into a huge sparse array requries as input of
%  each item its iy,ix,value triple and at the end, the size of this model
% see doc sparse
A=sparse([IE(:);IW(:);IN(:);IS(:);IT(:);IB(:)],...  % row    numbers of matrix coefficients
         [IW(:);IE(:);IS(:);IN(:);IB(:);IT(:)],...  % column numbers of matrix coefficients
        -[Cx(:);Cx(:);Cy(:);Cy(:);Cz(:);Cz(:)],...  % matrix coeffiicients themselves
         gr.Nxyz,...    % number of rows in system array
         gr.Nxyz,...    % number of columns in system array
         7*gr.Nxyz);    % total number of non-zero elements in system matrix

% Compute the diagonal values (but don't add them to the matrix yet,
% becauase we need the later on to compute the individal flow terms.
Adiag= -sum(A,2);

%% Mark which cells are active (Iact) which have to be computed (I), that is
%  the cells that are active but are not fixed-head cells, and, mark the cells
%  that are fixed-head cells, because we put them to the right hand side of the
%  matrix equation
IAct = find(IBOUND ~= 0);   % active cells
I    = find(IBOUND  > 0);   % active cells but not fixed head cells
Ifh  = find(IBOUND  < 0);   % active cells with fixed head

%% check if we must transpose to ensure vertical RHS arrays further down
mustTranspose = size(STRTHD(IAct),1)==1;

%% Let's go and solve

fprintf('Time steps:');

% We loop through time. We do this in backward fashion, just to prevent
% having to preallocate memory. This way it is done automatically. We'll
% return to the natural order when finished.
for i=Nt:-1:1
    fprintf('.'); % show that you'te busy.
    
    it = Nt-i+1;  % it = stress period nr, i is loop variable
    
    
    %% solve
    % this is the heart of everything. It says Ph = A\RHS with A the system
    % array and RHS the right-hand-side (with dimension L3/T].
    % spdiags put a diagonal into a sparse array. So here we put the
    % diagonal into our system array. But this diagonal is augmented by the
    % storag term. So storag terms + the already computed diagonal are
    % summed and inserted into the matrix A. But we only do this for the
    % active cells that are not fixed heads (I). We just ignore
    % the inactive cells and the fixed-head cells, they remain untouched.
    % On the right hand side are all the flow terms. First is the given
    % flows, FQ. Then we have the flows from the fixed heads we already
    % know. Finally, the known part of the storage term is added to the right
    % hand sied. Again, only the active cells that are not fixed-head cells
    % are involved. We thus get results for all active cells that are not fixed heads.
    % Solving is done using sparse matrix calculation, to save more than
    % 90% of memory requirements. Witout sparse matrix computations, the
    % maximum model size likely is 1000 times smaller.
    % So  this line, essentiall H = A\RHS is all that is needed to get the heads
    
    % In steady state Cs==0, so this is ok.
    if mustTranspose
        FQ=FQ'; Cs=Cs'; STRTHD=STRTHD';
    end
    H(i).values    = STRTHD;
    H(i).values(I) = spdiags(  Adiag(I)+Cs(I)/dt(it)/theta,   0,   A(I,I) ) \ ...
            (FQ(I)  -  A(I,Ifh)*STRTHD(Ifh)  +  Cs(I).*STRTHD(I)/dt(it)/theta);
                
    %% Next step is to get the different cel by cell flow terms and put those into
    % the budget struct under their correct label for later retrieval.
    
    if mustTranspose
        FQ=FQ'; Cs=Cs'; STRTHD=STRTHD';
        H(i).values = H(i).values';
    end
    
    % extend to end of time step and store (in case the computation is not fully implicit)
    H(i).values(IAct)=H(i).values(IAct)/theta-(1-theta)/theta*STRTHD(IAct);
    
    % reset STRTHD
    STRTHD = H(i).values;
    
    % add info to structs, compatibly with MODFLOW
    H(i).totim = t(it);     H(i).pertim = dt(it);
    H(i).time  = t(it);
    H(i).period= it;        H(i).tstep  = 1;
    
    if nargout>1        
        % step 1, store the labels
        B(i).label = labels;

        % preallocate to ensure correct size
        B(i).term{iCONST}  = gr.const(0);
        B(i).term{iRIGHT}  = gr.const(0);
        B(i).term{iFRONT}  = gr.const(0);
        B(i).term{iLOWER}  = gr.const(0);
        B(i).term{iWELLS}  = gr.const(0);
        B(i).term{iSTOR}   = gr.const(0);

        % step 2, compute internal balance of cell directly from the matrix and its
        % diagonal for all active cells, including the constant-head cells.
        B(i).term{iCONST}(IAct)  = spdiags(Adiag(IAct),0,A(IAct,IAct))*H(i).values(IAct);

        % we get the flow rightface, fontface and lowerface directly from the
        % heads and the intercell conductances.
        % We follow the modflow(1988) manual, p5.83-84 where the flow outward from
        % the cell across the mentioned faces is taken positive: (this gives
        % all the minus signs below:
        B(i).term{iRIGHT}(:,1:end-1,:)  = -Cx.*diff(H(i).values,1,2);
        B(i).term{iFRONT}(1:end-1,:,:)  = -Cy.*diff(H(i).values,1,1);
        B(i).term{iLOWER}(:,:,1:end-1)  = -Cz.*diff(H(i).values,1,3);

        % The rate of fixed head cells in this simple model is
        B(i).term{iWELLS}(IAct)  = FQ(IAct);

        % with the heads at the beginning and end of the time step known, we
        % can compute the storage rate
        B(i).term{iSTOR}(IAct) =-Cs(IAct).*(H(i).values(IAct)-STRTHD(IAct))/dt(it);   % Sstorage in time step m3 for cell

        % The real flow through the const head cells is then
        B(i).term{iCONST} = B(i).term{iCONST} - B(i).term{iWELLS} - B(i).term{iSTOR};
        % this should now be zero everywhere except at the constant head cells.
        % this zero is a good check of the water balance.
        B(i).budCheck = sum(B(i).term{iCONST}(I));

        B(i).time  = t(it);     B(i).pertim = dt(it);
        B(i).period= it;        B(i).tstp = 1;
    end
end
% flip back in order of increasing time
H = H(end:-1:1);

if nargout>1
    B = B(end:-1:1);
end

fprintf('%d time steps done\n',Nt);

if ~isempty(varargin)
    msgId = sprintf('mfLab:%s:vararginNotCompletelyUsed',mfilename);
    warning('on',msgId);
    warning(msgId,'%s: varargin not completely used',mfilename);
end

end

%% Checking size of input arrays
function check(gr,varargin)

    propNames = {'IBOUND','kx','ky','kz','Ss','STRTHD','FQ'};

    for i=1:numel(propNames)
        [Prop,varargin] = getProp(varargin,propNames{i},[]);
        propSize = size(Prop); propSize(end+1:3)=1;
        if  ~all(propSize == gr.size)
            error(['%s: %s must have size of grid, i.e.  [%d %d %d], not [%d %d %d].\n',...
                'REMEDY: use gr.const(%s) to generate a full array of constants.'],...
                mfilename,propNames{i},gr.size,propSize,propNames{i});
        end
    end

end

%%
function [H,B]=selfTest(n)
%SELFTEST tests function mFlow
%
% The selftest can be run by typing mFLow without arguments.
%
% TO 120422

if nargin==0 || isempty(n)
    for i=1:9
        [H,B] = mFlow(i);
    end
    return;
end

%% Inform user that selfTest is run
fprintf('\n');
msgId = sprintf('mfLab:%s:selfTest',mfilename);
warning('on',msgId);
warning(msgId,'%s Running selftest option %d',mfilename,n);
warning('off',msgId);


if n<1 || n>9
    error('%s: case number for self test must be between 1 and 9',mfilename);
end


%% Variables used by all examples below
k     = 25;    % [m/d] hydraulic conductivity
S     = 1e-3;  % [ - ] storage coefficient
c     = 360;   % [ d ] hydraulic resistance of top layer (if applicable)
d     = 10;    % [ m ] thickness of top layer (when applicable)
D     = 50;    % [ m ] thickness of aquifer
kD    = k*D;   % [L2/T] transmissivity of aquifer
lambda= sqrt(kD*c); % [ m ] spreading length of semi confined aquifer
Ss    = S/D;   % [1/m] specifi storage coefficient of aquifer
iWell = 2;     %  zone number in IBOUND indicating position of well screen
Qwell = -1200; % [m3/d] well extraction
xWell =     0; % [ m ] well coordinate in grid
yWell =     0; % [ m ] well coordinate in grid
t     = logspace(-4,2,61); % [d] simulation time (excluding initial time = 0)
N     = 0.01;  % [m/d];
iLay  =    1;  % [ - ] layer number for plotting in some cases

switch n
    case 1
        fprintf(2,'mFlow selfTest case %d: 3D De Glee flow model\n',n);
        iWell = 2;  % zone number for fixe flows
        Qwell = -1200;
        xWell = 0;
        yWell = 0;
        zWell = [-10 -30];

        %% 1: choose a grid
        xGr = logspace(0,4,30); xGr=[-xGr xGr];
        yGr = xGr;
        zGr = [0 -10 -50];

        gr = gridObj(xGr,yGr,zGr);

        %% set up IBOUND array and put zone number of well into it
        IBOUND        = gr.const(1);   % all cells active
        IBOUND(:,:,1) = -1;            % top layer fixed heads
        % put well at correct location into IBOUNND  with zoneNr == iWell
        IBOUND(hit(gr.xGr,xWell),  hit(gr.yGr,yWell),  fallsIn(gr.zGr,zWell)) = iWell;

        %% 3 choose time (if transient more than one time)
        t  =  logspace(-2,2,40);

        %% Set conductivities

        % In this case the vertical resistance of the top layer is given
        c   = 1000; % days resistance of top layer

        %% Allocate conductivities
        kx  = gr.const([NaN; 25]);

        % Replace top layer values with that computed from layer thickness and vertical hydraulic resistance
        kx(:,:,1) = gr.DZlay(:,:,1)/c/2;  % the 2 is because water travels through half DZtop only

        % use same values for ky and kz assuming:
        ky  = kx;  % isotropy horizontal
        kz  = kx;  % isotropy vertically

        Ss = gr.const(1e-3); Ss(:,:,1)=0;
        %% Specify initial heads to be all zeros
        STRTHD = gr.const(0);

        %% Specify fixed flows to be all zeros
        FQ = gr.const(0);
        % Make sure well is correct
        IW = IBOUND==iWell;  % use logical index

        % and compute FQ at well face accoring to the transmissivity of the
        % penetrated layers.
        FQ(IW) = Qwell * gr.DZlay(IW).*kx(IW)/sum(gr.DZ(IW).*kz(IW));

        %% Simulate and get output
        [H,B] = mFlow(gr,t,IBOUND,kx,ky,kz,Ss,STRTHD,FQ);

        %%
        hrange = ContourRange(H,50);
        iLay = 2;
        figure; hold on; xlabel('x [m]'); ylabel('y [m]');
        title(sprintf('mFlow selfTest case %d: %s selftest, head in layer %d',n,mfilename,iLay));
        contourf(gr.xc,gr.yc,H(end).values(:,:,iLay),hrange,'edgecolor','none');

        figure; hold on; xlabel('x [m]'); ylabel('drawdown [m]');
        title(sprintf('mFlow selfTest case %d: %s selftest, head different times',n,mfilename));
        iy = hit(gr.yGr,yWell);


        Ir = gr.xm>0;
        rm = gr.xm(Ir);

        D = gr.DZ(1,1,2);
        S = Ss(1,1,2)*D;
        kD= kx(1,1,2)*D;

        u      = (rm.^2).*S./(4.*kD);  % u without t
        lambda = sqrt(kD*c);
        rho    = rm/lambda;

        for it=1:numel(H)
            plot(log10(rm),H(it).values(iy,Ir,iLay));
            plot(log10(rm),Qwell/(4*pi*kD)*hantush(u/t(it),rho),'rx');
        end

        fprintf('%s\n','wait');

    case 2
        fprintf(2,'mFlow selfTest case %d: Example Hantush 3D (of quasi 3D)\n',n);
        
        % generate a suitable grid
        xGr = logspace(-1,4,60); xGr =[ -xGr xGr]; yGr = xGr; zGr = [d 0 -D];  gr = gridObj(xGr,yGr,zGr);

        % expand the sclar variables to the size of the grid
        K = gr.const(k); SS=gr.const(Ss); STRTHD=gr.const(0); FQ=gr.const(0); IBOUND=gr.const(1);

        K( :,:,1) = d/c/2;  % compute vertical conductivity of top layer form resistance and thickness
        SS(:,:,1) = 0;      % set specific storage of top layer to zero, to match analytical solutio

        iLay = 2;           % layer with well screen
        IBOUND(:,:,1) = -1; % top layer gets fixed head
                            % put well in the grid by setting zone number to mark
                            % its position in the grid
        IBOUND(  hit(gr.xGr,xWell), hit(gr.yGr,yWell), iLay) = iWell;

        FQ(IBOUND==iWell) = Qwell; % insert the well

        % Run the model and get its output
        [H,B] = mFlow(gr,t,IBOUND,K,K,K,SS,STRTHD,FQ);

        % Show the results along a radial line for which we select the x-axis for xm>0
        IR = gr.xm>0;   % logical indices suitable a radius
        rm = gr.xm(IR); % radius rm

        u      = rm.^2*S/(4*kD); % u without t
        rho    = rm/lambda; 

        figure; hold on;
        title(sprintf('mFlow selfTest case %d: example Hantush',n));
        xlabel('r [m]'); ylabel('head [m]');

        iy = hit(gr.yGr,yWell);  % row through well
        for it=1:2:numel(H)
            plot(rm,H(it).values(iy,IR,iLay),'b'); % numerical
            % compute analytic values
            for ir=numel(rho):-1:1
                Htsh(ir) = hantush(u(ir)/t(it),rho(ir));
            end
            plot(rm,Qwell/(4*pi*kD)*Htsh,'rx'); % analytical
        end
        set(gca,'xscale','log');

        zonebudget(B);
    case 3
        fprintf(2,'mFlow selfTest case %d: Hantush axially symmetric\n',n);
        
        % Here we will run the model in axially symmetric mode. We only need a
        % cross section in this case. We set yGr = [], which is interpreted as 
        % [-0.5 0.5], i.e. 1 m wide. And we set AXIAL=true when generating the
        % grid:

        rGr = logspace(-1,4,60); yGr = []; zGr = [d 0 -D];  gr = gridObj(rGr,yGr,zGr,'AXIAL',true);

        % As in the previous example:
        K = gr.const(k); SS=gr.const(Ss); STRTHD=gr.const(0); FQ=gr.const(0); IBOUND=gr.const(1);
        K( :,:,1) = d/c/2;
        SS(:,:,1) = 0;

        iLay = 2; % layer with well screen
        IBOUND(:,:,1) = -1;
        IBOUND(1,1,iLay) = iWell;
        FQ(IBOUND==iWell) = Qwell;

        % Run the model and get its output
        [H,B] = mFlow(gr,t,IBOUND,K,K,K,SS,STRTHD,FQ);

        % AXIAL=true, so we can use gr.rm directly
        u      = gr.rm.^2*S/(4*kD); % u without t
        rho    = gr.rm/lambda; 

        figure; hold on;
        title(sprintf('mFlow selfTest case %d: Example Hantush',n));
        xlabel('r [m]'); ylabel('head [m]');
        
        iy=1;
        for it=1:2:numel(H)
            plot(gr.rm,H(it).values(iy,:,iLay),'b'); % numerical
            % compute Hantush
            for ir=numel(rho):-1:1
                Htsh(ir) = hantush(u(ir)/t(it),rho(ir));
            end
            plot(gr.rm,Qwell/(4*pi*kD)*Htsh,'rx');   % analytical
        end
        set(gca,'xscale','log');

        zonebudget(B);
        
    case 4
        fprintf(2,'mFlow selfTest case %d: Theis, a well in an infinite aquifer\n',n);
        
        % Grid as before, but we need only one layer in this case

        xGr = logspace(-1,5,60); xGr =[ -xGr xGr]; yGr = xGr; zGr = [0 -D];  gr = gridObj(xGr,yGr,zGr);

        K = gr.const(k); SS=gr.const(Ss); STRTHD=gr.const(0); FQ=gr.const(0); IBOUND=gr.const(1);

        iLay=1; % well is now in layer 1
        IBOUND( hit(gr.xGr,xWell),hit(gr.yGr,yWell),iLay) = iWell;
        FQ(IBOUND==iWell) = Qwell;

        [H,B] = mFlow(gr,t,IBOUND,K,K,K,SS,STRTHD,FQ);

        IR = gr.xm>0;   % logical indices suitable a radius
        rm = gr.xm(IR); % radius

        u = rm.^2*S/(4*kD); % u without t

        figure; hold on; xlabel('r [m]'); ylabel('head [m]');
        title(sprintf('mFlow selfTest case %d: Example Theis',n));

        iy = hit(gr.yGr,yWell);  % row through well
        for it=5:5:numel(H)
            plot(rm,H(it).values(iy,IR,iLay),'b');              % numeric
            plot(rm,Qwell/(4*pi*kD)*expint(u/t(it)),'rx');   % analytic
        end
        set(gca,'xScale','log');

        zonebudget(B);
        
    case 5
        fprintf(2,'mFlow selfTest case %d: Theis axially symmetric\n',n);
        
        % Again the grid has one layer and AXIAL=true is set when generating the grid

        rGr = logspace(-1,4,60); yGr = []; zGr = [0 -D];  gr = gridObj(rGr,yGr,zGr,'AXIAL',true);

        K = gr.const(k); SS=gr.const(Ss); STRTHD=gr.const(0); FQ=gr.const(0); IBOUND=gr.const(1);

        iLay = 1;
        IBOUND(1,1,iLay) = iWell;
        FQ(IBOUND==iWell) = Qwell;

        % Run the model and get its output
        [H,B] = mFlow(gr,t,IBOUND,K,K,K,SS,STRTHD,FQ);

        u      = gr.rm.^2*S/(4*kD); % u without t

        figure; hold on;
        title(sprintf('mFlow selfTest case %d: Example Theis, axisymmetriC',n));
        xlabel('r [m]'); ylabel('head [m]');
        
        iy = 1;
        for it=1:2:numel(H)
            plot(gr.rm,H(it).values(iy,:,iLay),'b');
            plot(gr.rm,Qwell/(4*pi*kD)*expint(u/t(it)),'rx')
        end
        set(gca,'xscale','log');

        zonebudget(B);
        
    case 6
        fprintf(2,'mFlow selfTest case %d: Dupuit: steady state well with fixed boundary at r=R\n',n);

        xGr = logspace(-1,log10(3500),60); xGr =[ -xGr xGr]; yGr = xGr; zGr = [0 -D];  gr = gridObj(xGr,yGr,zGr);

        K = gr.const(k); SS=gr.const(Ss); STRTHD=gr.const(0); FQ=gr.const(0); IBOUND=gr.const(1);

        iLay = 1;
        R    = 1750;  % radius at which drawdown is zero
        IBOUND(gr.RM>=R)=-1;  % set fixed head at x,y>R
        IBOUND( hit(gr.xGr,xWell),hit(gr.yGr,yWell),iLay) = iWell;
        FQ(IBOUND==iWell) = Qwell;

        [H,B] = mFlow(gr,t(1),IBOUND,K,K,K,SS,STRTHD,FQ);

        figure; hold on;
        title(sprintf('mFlow selfTest case %d: Steady state Dupuit example',n));
        xlabel('r [m]'); ylabel('head [m]');
        
        iy = hit(gr.yGr,yWell);
        IR = gr.xm>0; rm=gr.xm(IR);
        plot(rm,H(1).values(iy,IR,iLay),'b');
        plot(rm,Qwell/(2*pi*kD)*log(R./rm),'ko');
        set(gca,'xScale','log');

        zonebudget(B);
    case 7
        fprintf(2,'mFlow selfTest case %d: Circular Island with Precipitation\n',n);

        xGr = -2000:25:2000; yGr = xGr; zGr = [0 -D];  gr = gridObj(xGr,yGr,zGr);

        K = gr.const(k); SS=gr.const(Ss); STRTHD=gr.const(0); FQ=gr.const(0); IBOUND=gr.const(1);

        R    = 1750;
        iLay = 1;
        IBOUND(gr.RM>=R)=-1;

        N = 0.001;
        FQ(:,:,1) = N*gr.AREA;

        % Run the model and get its output
        [H,B] = mFlow(gr,t(1),IBOUND,K,K,K,SS,STRTHD,FQ);

        figure; hold on;
        title(sprintf('mFlow selfTest case %d: Steady state Dupuit example',n));
        xlabel('r [m]'); ylabel('head [m]');
        
        iy = hit(gr.yGr,0);
        IR = gr.xm>0; rm = gr.xm(IR);
        plot(rm,H(1).values(iy,IR,iLay),'b');
        plot(rm,N/(4*kD) * (R^2-rm.^2),'rx');

        zonebudget(B);
        
    case 8
        fprintf(2,'mFlow selfTest case %d: Linear with recharge along x-axis\n',n);

        xGr =  [0.01 0:100:1000 1000.01];
        yGr =  1000:-100:0;
        zGr =  [0 -D];
        gr = gridObj(xGr,yGr,zGr);

        K     = gr.const(k);
        SS    = gr.const(Ss);
        STRTHD= gr.const(0);

        IBOUND= gr.const(1);
        IBOUND(:,end,:)=-1;

        q  = 1; % m2/d 
        FQ = gr.const(0);
        FQ(:,:,1) = N*gr.AREA;
        %FQ(:,1,:) = q * gr.dy;

        % Run the model and get its output
        [H,B] = mFlow(gr,t(1),IBOUND,K,K,K,SS,STRTHD,FQ);

        figure; hold on;
        title(sprintf('mFlow selfTest case %d: Steady state Dupuit example',n));
        xlabel('r [m]'); ylabel('head [m]');
        iy=1;
        plot(gr.xm,H(1).values(iy,:,iLay),'b');
        plot(gr.xm, q/kD*(sum(gr.dx)-gr.xm),'rx');
        plot(gr.xm, N/(2*kD)*(sum(gr.dx).^2-gr.xm.^2),'rx');

        zonebudget(B);
    case 9
        fprintf(2,'mFlow selfTest case %d: Linear with recharge along y-axis\n',n);

        xGr =  [0.01 0:100:1000 1000.01];
        yGr =  [0.01 0:100:1000 1000.01];
        zGr =  [0 -D];
        gr = gridObj(xGr,yGr,zGr);

        K     = gr.const(k);
        SS    = gr.const(Ss);
        STRTHD= gr.const(0);

        IBOUND= gr.const(1);
        IBOUND(1,:,:)=-1;

        q  = 1; % m2/d 
        FQ = gr.const(0);
        FQ(:,:,1) = N*gr.AREA;
        %FQ(:,1,:) = q * gr.dy;

        % Run the model and get its output
        [H,B] = mFlow(gr,t(1),IBOUND,K,K,K,SS,STRTHD,FQ);

        figure; hold on;
        title(sprintf('mFlow selfTest case %d: Steady state Dupuit example',n));
        xlabel('r [m]'); ylabel('head [m]');
        
        ix = 1;
        plot(gr.ym,H(1).values(:,ix,iLay),'b');
        plot(gr.ym, q/kD*(sum(gr.dy)-gr.ym),'rx');
        plot(gr.ym, N/(2*kD)*(sum(gr.dy).^2-gr.ym.^2),'rx');

        zonebudget(B);
        
        fprintf(2,'\nDone\n');
        
    otherwise
        fprintf(2,'Unknown option use 1..9\n');
end

end


    
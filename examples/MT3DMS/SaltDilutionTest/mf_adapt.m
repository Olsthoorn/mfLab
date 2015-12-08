% Salt dilution test simulated with MT3D

% TO 091029 091204
%
% We simulate the flow in a pipe with a square cross section which is
% connected to a box (cube) and a subsequent pipe with the same XSection.
% The right end of the model has a fixed head.
% The left end has a fixed flow boundary.
% Other boundaries are all closed.
% Initial concentrations are zero, the flow is steady state.
%
% To simulate flow in a pipe and box we set porosity to 1.
% The idea is to inject salt during a short period and measure the
% salinity curve as it passes by observation points in the system.
% We then use these curves to estimate the total flow through the system,
% which is the essence of a salt dilution test.
%
% This configuration resembles the flow in extraction galleries in the
% Amsterdam Water Supply Dunes in which out salt dilution tests were
% carried out to measure the flow along the 300 m long extraction gallery
% at its obsrvation points, which are 20 cm wide pipes that are connected
% to the gallery and extend to above ground surface.
% 
% Real flow in the gallery is, of course, turbulent, which we cannot simulate
% with a groundwater transport model. However, we may see whether putting the
% salinity sensors at different places in the box makes more or less sense.
% In this simulation the flow situation is completely known and therefore
% suitable for testing the numeric model, which is the purpose of this example.
%
% Three stress periods are simulated:
%   First 5 seconds with zero concentration input and then the salt is dosed
%   during 5 seconds, after which the input concentration is again zero.
%   After 5 seconds salt is dosed during 5 seconds.
%   Subsequently, the salt is washed out by the flow.
%
% Output is produced at every second.
%
% One may want to experiment with dispersion etc.
% Observation points may be added, removed or changed to explore the
% concentration curvers at other places.
% points.
%
% TO 091204
%
%% recharge on a rectangular area with fixed head boundary at left and right
clear variables;
close all;

basename='SaltTest';

BACKGROUND = true;
GREP = 'PERIOD';

%% The model grid
n    = 15; % number of cells used in cross section
NROW = n;
NCOL = 3*n;
NLAY = n;

ribLength = 0.04; % m ribLength of model cell
Rpipe     = 0.2;  % m radius of pipe that is connected to box

xGr=(0:NCOL)   *ribLength;
yGr=(0:NROW)   *ribLength;
zGr=(NLAY:-1:0)*ribLength;

gr=gridObj(xGr-mean(xGr),yGr-mean(yGr),zGr-mean(zGr));

%% fixed head boundaries. Are interpreted as local point heads in fdm2dens
STRTHD = gr.const(0); % initial heads
HK     = gr.const(1); % conductivity (immaterial here)
VK     = gr.const(1); % same, vertially

IBOUND = gr.const(0); % determine which model cells to copute, keep and exclude

R = sqrt(gr.YM.^2+gr.ZM.^2);

IBOUND(R<=Rpipe)       = 1;
%IBOUND(6:11,  :  ,6:11)= 1;    % activate the pipe cells
IBOUND( :  ,16:31, :  )= 1;    % activate the box  cells
IBOUND(6:11, end ,6:11)=-1;    % set right boundary to fixed heads

ICBUND = gr.const(1); % same as IBOUND but for concentrations
PEFF   = gr.const(1); % effective porosity
STCONC = gr.const(0); % initial concentration of tracer

%% Visualisation of cross section by plotting IBOUND itself

figure; title('Cross section of gallery with shaft in the center');
xlabel('xGr [m]'); ylabel('z [m]');
contourf(gr.xc,gr.zc,XS(IBOUND(8,:,:)),[0.5,-0.5]);  % show layers

%% Get NPER from workbook
[PERparnams,PERparvals,NPER]=getPeriods(basename);
% To specify this injection concentration in a most general fashion,
% we add a new variable to the PER sheet, called cInj and specify the
% tracer concentrations for all stress periods (either zero or the value we
% desire. Having done so, we can retrieve the injected tracer concentration
% as follows:
CInj=PERparvals(:,strmatchi('CInj',PERparnams));

% The duration of the stress periods comes from
PERLEN=PERparvals(:,strmatchi('PERLEN',PERparnams));

% Ok but how much mass have we injected?
% To compute this, we use the water velocity, the cross section, which is
% given here:
v0 =0.15;       % m/s   water velocity

% To compute the cross section of the pipe, we make use of IBOUND, which
% has value 1 for active and 0 for ineactive cells. There are different
% ways to compute this cross section
% One way is to turn a cross seciton toward us (align z vertically and rows
% from left to right while taking a cross section (a given column, for
% instance the second one). This is done using permute, which yield the
% cross section plane wigh zeros for the inactive cells and ones for the
% active cells. To show it try:
figure
spy(YS(IBOUND(:,3,:))); xlabel('iy'); ylabel('iz');

% If we multiply this element by element (.*) with the cross section area
% of the individual cells in the plane, obtained by the matrix
% multiplicaton (*) of DELZ vertically and DELY horizontally, wer're done:
A=sum(sum(YS(IBOUND(:,3,:)).*(gr.dz(:)*gr.dy')));
% this is somewhat complicated. There is a somewhat easier way to compute
% the cross section area as is shown in the next section

% Ass all cInj in stress periods without tracer are zero, we can just
% compute the entire injected mass over all stress periods as follows:
M=v0*A*sum(CInj.*PERLEN);  % g if conc is in g/cm3 here 

%% Mass loading at single entry cell or injection with a given

% We may inject the tracer by means of a uniform concentration over the entire cross
% section. Or we may inject pure tracer at one or more of the entry cells by setting
% the concentrataion at these points (ITYPE=2 (WEL) in SSM) or by injection mass
% at these points (ITYPE=15) in the SSM module of MT3DMS.
%
% Both ITYPE are simulated by means of the switch massLoading in this mf_adapt m-file:
% If on  (true), then mass is injected at a single cell at the entry face of the model,
% if off (false) a constant tracer concentration across the entry face is maintaed.
%
% The same parameters are used in both cases. In the case of massLoading CInj
% is interperted as g/s instead of concentration
%

massLoading = true;   % use false for "off" and true for "on"

%% concentration during the second stress period across the entire zone

% Which are the left hand side cell in the pipe??
iLeft = 7; % zone number for left side of model within pipe
IBOUND(IBOUND>0 & gr.XM<gr.xGr(2)) = iLeft;

% Alternative way to compute the injection mass
I = find(IBOUND==iLeft);
A = sum(gr.DY(I).*gr.DZ(I));

% Total discharge of model given v0 and configuration
Q0=v0*gr.DY(I).*gr.DZ(I); % m3/s

switch massLoading    % we select a random cell from the model entry face
    case true
        %% Choose some random cells at the model entry face as injection points.    
        irand= round(rand(3,1)*numel(I-1))+1;
        % Compute the tracer entry mass in each stress period
        M = sum(CInj.*PERLEN); % where now CInj is interpreted as mass inje 
        [WEL,PNTSRC] = bcnZone(basename,'MLS',IBOUND,{iLeft Q0},{'cInj'});
    case false
        % If no massLoading, treat CInj as uniform conc over the entry face
        M = v0*A*sum(CInj.*PERLEN);
        [WEL,PNTSRC] = bcnZone(basename,'WEL',IBOUND,{iLeft Q0},{'cInj'});
end

%[MLL,PNTSRC] = bcnZone(basename,'MML',IBOUND,{iLeft Q},{'cInj'});

save underneath Q0 M v0 iLeft

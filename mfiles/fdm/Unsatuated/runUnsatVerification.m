% Run unsat1

close all
clear

%% Determine what to show
showMeteo = true;
showOut   = true;
showBudget= true;
uitzakking= false;

sheetNm = 'params';

%% get Meteo
meteoDir = '/Users/Theo/GRWMODELS/mflab/examples/TimeSeriesAnalysis/KNMI-data-files/';
load([meteoDir 'TPE_Rotterdam344']);  % loads TPS for Rotterdam

%% select simulation period

verificationItem = 10;

switch verificationItem
    case 1 % verify fixed head boundary condition
        % set fixeH in workbook to some value not NA()
        % seepage to zero
        P = 0.000 * ones(90,1) * [0 1 0 2 0];
        E = 0.000 * ones(90,1) * [0 0 1 0 2];
        TPE = [datenum(2015,1,1)+cumsum(ones(size(P(:)))) P(:) E(:)];
        % verify head and flow
    case 2 % Constant infiltration
        % Run with seepage=0
        % FH = -7.9 m at zm=-7.94 m
        % Initially, h=0 all along the profile.
        % Show that h=0 at all model cells
        P = 0.01298 * ones(90,1) * [1 1 1 1];
        P = P(:);
        E = zeros(size(P));
        TPE = [datenum(2015,1,1)+cumsum(ones(size(P(:)))) P(:) E(:)];
    case 3 % Constant extraction through top of model
        % Same situation as before, but using negative P to extract water
        % from the top of the model. Make sure maxIC=0 (no interception)
        P = -0.001 * ones(90,1) * [0 1 1 1 1 1 1];
        E =  0.000 * ones(90,1) * [0 0.5 1 2 4 8 16];
        TPE = [datenum(2015,1,1)+cumsum(ones(size(P(:)))) P(:) E(:)];       
    case 4 % constant evaporation
        % Same as before, but now extraction is from the rootzone. Also the
        % E was raised in steps from 1.0 to 4.5 mm/d, causing E-reduction
        % to become active through the Feddes(1978) reduction.
        P =  0.000 * ones(90,1) * [0 1 1 1 1 1 1];
        E =  0.0010 * ones(90,1) * [0 1 1 1 1 1 1];
        TPE = [datenum(2015,1,1)+cumsum(ones(size(P(:)))) P(:) E(:)];
        
        % This situation with E=1.0 mm/d was also used to verify that
        % the total inflows match the storage change, QstoSS, QstoTheta
        % The methods to compute the storage were compared.
    case 5 % verify drain at model bottom
        % Zet fixedH = NA();
        P =  0.005 * ones(90,1) * [0 1 1 1 1 1 1];
        E =  0.000 * ones(90,1) * [0 1 1 1 1 1 1];
        TPE = [datenum(2015,1,1)+cumsum(ones(size(P(:)))) P(:) E(:)];
        
    case 6 % verify surface runoff by drain at top
        P =  0.05 * ones(90,1) * [0 1 2 3 4 5 6];
        E =  0.000 * ones(90,1) * [0 1 1 1 1 1 1];
        TPE = [datenum(2015,1,1)+cumsum(ones(size(P(:)))) P(:) E(:)];
       
    case 7 % Verify ditch
        % ground surface -4 m, hDitch -4.5 m, no fixed head, no seepage
        % wEx = 2 d, wIn = 10 d.
        P =  0.000 * ones(90,1) * [0 1 1 1 1 1 1];
        E =  0.001 * ones(90,1) * [0 1 1 1 1 1 1];
        TPE = [datenum(2015,1,1)+cumsum(ones(size(P(:)))) P(:) E(:)];

        
    case 8 % verify seepage face flow above dich
        % show that the downward progression of the wetting front matches
        % the analyical solution by Greem Ampt.
        P = 0.001 * ones(90,1) * [0 10 0];
        E = zeros(size(P));
        TPE = [datenum(2015,1,1)+cumsum(ones(size(P(:)))) P(:) E(:)];       
    case 9 % sharp front infiltration
        sheetNm = 'paramsDiep';
        % with a closed model show that decline of head matches the
        % evaportranspiration and the storage coefficient.
        % also show the relation between Ss and delta theta
        P = 0.001 * ones(90,1) * [1 2 4 8 16 32 64 128 256 512];
        E = 0.000 * ones(90,1) * [0 1 1 1  1  1  1   1   1   1];
        TPE = [datenum(2015,1,1)+cumsum(ones(size(P(:)))) P(:) E(:)];  
    case 10 % sharp front infiltration
        sheetNm = 'paramsDiep';
        % with a closed model show that decline of head matches the
        % evaportranspiration and the storage coefficient.
        % also show the relation between Ss and delta theta
        P = 0.050 * ones(5,1) * [0 1 zeros(1,70) [0 1 zeros(1,70)] [0 1 zeros(1,70)]];
        E = 0.004 * ones(5,1) * [1 1 ones( 1,70) [1 1 ones( 1,70)] [1 1 ones(1,70)]];
        TPE = [datenum(2015,1,1)+cumsum(ones(size(P(:)))) P(:) E(:)];  
    otherwise
        error('Unknown verificationItem');
end


%TPE(:,3) = 0;

%% meteo Rotterdam
if showMeteo==true
    figure; axes('nextPlot','add','xGrid','on','yGrid','on','fontSize',12);
    xlabel('2010'); ylabel('m/d'); title('P and E in Rotterdam, station 344');

    
    bar(TPE(:,1),TPE(:,2),'b');
    plot(TPE(:,1),TPE(:,3),'r');
    datetick;
    
    % cumulative meteo Rotterdam vanaf 1 april 2010
    I = TPE(:,1)>=datenum(2010,4,1);
    figure; plot(TPE(I,1),cumsum(TPE(I,2)),'b',...
            TPE(I,1),cumsum(TPE(I,3)),'r',...
            TPE(I,1),cumsum(TPE(I,3)-TPE(I,2)),'k');
        datetick;
    xlabel('2010'); ylabel('m'); title('P and E en E-P cumulatief in Rotterdam, station 344');
    legend('P','E','E-P');
    grid on;
end


runVerificationCase

showUnsatResults

unsatBudget

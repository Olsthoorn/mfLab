%% Analyzing output of the model
% TO 190727

if exist('AFTERMFSETUP','var'), PEST=1; else PEST=0; end

if ~PEST
    load('name.mat') % loads name.mat, which only contains the variable "basename"
    load(basename);  % having retrieved the baasename value, we can load
end

load underneath

H=readDat([basename,'','.hds']);  % use readDAT to read the heads  output file
H=maskHC(H,1000);

hObs=H(1).values(IObs);

if PEST
%% Printing the model output file for PEST
    fid=fopen([basename,'.out'],'w');
    fprintf(fid,'%15.5g\n',hObs');
    fclose(fid);
    
    exit;  % close matlab
end

if ~PEST
    NLAY=size(IBOUND,3); % get the size of the model

    [xGr,yGr,xm,ym,DELX,DELY,NCOL,NROW]=modelsize(xGr,yGr); % xGr and yGr were saved in mf_adapt

    figure
    for iLay=1:NLAY    % plot 3D surface of every layer of model
        subplot(NLAY,1,iLay)
        h=surf(xm/1000,ym/1000,H.values(:,:,iLay)); hold on;
        title(sprintf('Example 1 of mf2k Open file Report 00-92 p89 [Layer %d]',iLay));
        xlabel('x [1000ft]');  ylabel('y [1000ft]');
    end

    % This shows that the difference between te layers is small
    %% Contouring heads

    figure
    iLay=3;
    [c,hdl]=contour(xm,ym,H.values(:,:,iLay));  % default contouring
    clabel(c,hdl);                          % default labels
    title(sprintf('Head contours, Open file Report 00-92 p89 [Layer %d]',iLay));
    xlabel('x [1000ft]');  ylabel('y [1000ft]');

    %% Specific contouring
    phirange=ContourRange(H,25);

    figure                               % new picture
    for iLay=1:NLAY
        subplot(1,NLAY,iLay);
        [c,hdl]=contourf(xm/1000,ym/1000,H.values(:,:,iLay),phirange);
        clabel(c,hdl);
        axis('equal'); % make both axes same scale
        axis('tight'); % no white area around axes
        caxis(phirange([1 end]));  % fix color scale to these values
        colorbar;      % place colorbar but now with known scale
        xlabel('x [1000ft]'); if iLay==1, ylabel('y [1000ft]'); end
        title(sprintf('period=%d, tstp=%d, layer=%d',...
                       H(end).period,H(end).tstp,iLay));      
    end

    %% Read unformatted budget file and mask noflow cells if they exist

    B=readBud([basename,'','.bgt']);  B=maskHC(B,IBOUND);

    totDRN=sum( B(1).term{ strmatch('DRAINS',      B(1).label) }(:) );
    totRCH=sum( B(1).term{ strmatch('RECHARGE',    B(1).label) }(:) );
    totCHD=sum( B(1).term{ strmatch('CONSTANTHEAD',B(1).label) }(:) );
    totWEL=sum( B(1).term{ strmatch('WELLS'       ,B(1).label) }(:) );

    fprintf(['total recharge is             = %12.2f ft^3/d\n',...
             'total nwt well discharge      = %12.2f ft^3/d\n',...    
             'total constant head discharge = %12.2f ft^3/d\n',...
             'total drain discharge         = %12.2f ft^3/d\n',...
             '                                ------------------\n',...
             'overall total (balance)       = %12.f  ft^3/d\n'],...
                totRCH,totWEL,totCHD,totDRN, totRCH+totWEL+totCHD+totDRN);

    %% In and outflow through constant head cells
    CH=B.term{strmatch('CONSTANTHEAD',B.label)};
    QCHin =sum(CH(CH(:)>0));
    QCHout=sum(CH(CH(:)<0));

    fprintf(['Total inflow  across constant head cells    =  %12.2f ft^3/d\n',...
             'Total outflow across constant head cells    =  %12.2f ft^3/d\n',...
             '                                               -----------------\n',...
             'Total net inflow across constant head cells =  %12.2f ft^3/d\n'],...
             QCHin,QCHout,QCHin+QCHout);

    %% Zonebudget
    iz=reshape(1:3,[1,1,3]);  % one zone per layer
    Zonearray=iz(ones(NROW,1),ones(1,NCOL),:);
    zonebudget(B,Zonearray,1);
    zonebudget(B,Zonearray,2);
    zonebudget(B,Zonearray,3);
    zonebudget(B,Zonearray,[1 3]);
    zonebudget(B,Zonearray,[1 2 3]);
end

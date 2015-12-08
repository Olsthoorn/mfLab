function writeMPTH6(~,pth)
%WRITEMPTH6 writes MOPATH6's simulation input file
%
% Example:
%    writeMPTH(basename,pth)
%
% See also: writeMPBAS
%
% TO 130121

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

GRID=1;  % currently, MODPATH allows only one grid (Local Refinement not implemented in MODPATH6 (2012))

ctrlRecordDesired = true;
freeFormat        = true;

fid=fopen(pth.simulationFile,'wt');

%0.
fprintf(fid,'# Matlab  MODPATH v6.0 simulation file %s %s\n',mfilename,datestr(now));
fprintf(    '# MODFLOW MODPATH v6.0 simulation file %s %s\n',mfilename,datestr(now));

%1.
fprintf(fid,'%s\n',pth.nameFile);

%2
fprintf(fid,'%s\n',pth.listingFile);

%3
fprintf(fid,repmat(' %d',[1,12]),...
    pth.simulationType,...            % 1 endpoint, 2 pathline, 3 is timeseries simulation
    pth.trackingDirection,...         % 1 forward, 2 backward
    pth.weakSinkOption,...            % always use 1, i.e. allow passing weak sinks
    pth.weakSourceOption,...          % always use 1, i.e. allow passing weak sources
    pth.referenceTimeOption,...       % always use 1 i.e. specify a reference time
    pth.stopOption,...                % always use 1 stop at end or start of Modflow Simulation
    pth.particleGenerationOption,...  % always use 2 (read from startingLocationsFile)
    pth.timePointOption,...           % we will always specify time points or zero
    pth.budgetOutputOption,...        % always use 1 (no budget checking)
    pth.zoneArrayOption,...           % always use 1, no zone data is read
    pth.retardationOption,...         % always use 2, i.e. if present
    pth.advectiveObservationsOption); % always use 1, because it's nonsense.

    fprintf(fid,'     simTyp trckDir wkSnk strSnk refTimeOp stopOp ptclGenOp tPntOp bgtOutOp zoneOp retardOp advObsOp\n');

    %4 endpoint file  (implied by mfLab)
    fprintf(fid,'%s\n',pth.endPointsFile);
    
    switch pth.simulationType
        case 1
        case 2
            %5
            fprintf(fid,'%s\n',pth.pathLinesFile);
        case 3
            %6
            fprintf(fid,'%s\n',pth.timeSeriesFile);
            if pth.advectiveObservationsOption == 2
                %7
                fprintf(fid,'%s\n',pth.advectiveObservationsFile);
            end
    end
    
    switch pth.referenceTimeOption
        case 1
            %8
            fprintf(fid,'%g     referenceTime\n',pth.referenceTime);
        case 2
            %9
            fprintf(fid,'%d  %d  %g     Period Step TimeFraction\n',...
                pth.period,pth.step,pth.timeFraction);
    end
    
    %10
    switch pth.stopOption
        case 1
        case 2
        case 3
            fprintf(fid,'%g     stopTime\n',pth.stopTime);
    end
    
    %% Particle starting locations
    switch pth.particleGenerationOption
        case 1, % write the particle locations file here
            pth.particleGroup.writeStartingLocations(fid);
        case 2, % use a separate file for the particle generation
            
            fprintf(fid,'%s\n',pth.startingLocationsFile);
            
            pth.particleGroup.writeStartingLocations(pth.startingLocationsFile);
            %22
        otherwise
            error('particleGenerationOption must be 1..2 see worksheet MPATH');
    end
    
    %% timePoints
    switch pth.timePointOption
        case 1
        case 2
           %23
           fprintf(fid,'%d     # of timePoints\n',pth.timePointCount);
           %24
           fprintf(fid,'%g     releaseTimeIncrement\n',pth.releaseTimeIncrement);
        case 3
           %23
           fprintf(fid,'%d     # of timePoints\n',pth.timePointCount);
           %25
           fprintf(fid,' %g',pth.timePoints(1:pth.timePointCount));
           fprintf(fid,'\n');
        otherwise
            error('timePointOptions must be 1..3, see worksheet MPATH');
    end

    
    switch pth.budgetOutputOption
        case {1,2}  %do nothing
        case 3
            %26
            fprintf(fid,'%d     # of BudgetCells\n',numel(pth.BUDGETCELLS,1));
            %27
            fprintf(fid,sprint('%d %%4d %%4d %%4d',GRID),...
                cellIndices(pth.BUDGETCELLS,pth.dims,'LRC')');                
        case 4
            %28
            fprintf(fid,'%s\n',mpbas.traceFile);
            %29
            fprintf(fid,'%d\n',mpbas.traceId);
    end
    

    switch pth.zoneArrayOption
        case 1
        case 2
            %30
            fprintf(fid,'%d      stopZone\n',pth.stopZone);
            %31
            for iLay=1:numel(pth.LAYCBD)
                warray(fid,pth.ZONE(:,:,iLay),0,'(25I3)',sprintf('ZONE{%d}',iLay),...
                    ctrlRecordDesired,freeFormat);
            end
    end    
    
    switch pth.retardationOption
        case 1
        case 2
            %32
            icbd=0;
            for iLay=1:gr.NLay
                warray(fid,pth.RETARDATION(:,:,iLay),0,'(10E15.6)',sprintf('RetardationFactor{%d}',iLay),...
                    ctrlRecordDesired,freeFormat);
                if pth.LAYCBD(iLay)>0
                    icbd=icbd+1;
                    %33 [RetardationFactorCB(NCOL,NROW)] if any GRID.LAYCBD>0
                    warray(fid,pth.RETARDATIONCB(:,:,icbd),0,'(10E15.6)',sprintf('RETARDATIONCB {%d}',icbd),...
                        ctrlRecordDesired,freeFormat);
                end
            end
    end
    
fclose(fid);

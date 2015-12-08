function LOCATIONS = startLoc(o,zoneArray,zoneVals)
    % gridObj/startLoc: startloc = gr.startLoc(basename,zoneArray,zoneVals)
    %
    % LOCATIONS are starting locations for MODPATH
    % zoneArray is an aray with one line per zone having the following items
    % { zoneNr nLoc iFace t_released }
    % nLoc determines number of points to be distributed over iFace. nLoc can be
    % [Nx Ny Nz] or just n which means Nx=Ny=Nz=n.
    % iFace determines the cell face at which the particles are released
    %  iFace= 1=West 2=East, 3=Front, 4=Back, 5=Bottom, 6=Top 0=Interior distribution.
    % iFace      multiple iFace may be used to set them simultaniously
    % 
    % t_released may be defined per zone
    %
    % LOCATIONS in the workspace will be picked up by MODPATH when invoked
    %
    % TO 120528
    
    nZones = size(zoneVals,1);
    
    LOCATIONS{nZones,1} =[];
    
    for izone=1:nZones
        iZone = zoneVals{izone,1};
        n     = zoneVals{izone,2}; % [Nx Ny Nz]
        iFace = zoneVals{izone,3};
        
        if size(zoneVals,2)<4
            t_released = 0;
        else
            t_released = zoneVals{izone,4};
        end
        
        CRL = cellIndices(zoneArray==iZone,o.size,'CRL');
        
        relLoc = distribute(n,iFace);
        
        nCells  = size(CRL,1);        % Nr of cells in izone
        nPoints = size(relLoc,1);   % Nr of starting points per cell
        
        % allocate for starting points of this zone
        LOCATIONS{izone,1}=NaN(nCells*nPoints,10);
        
        % get the starting points
        k=0;
        for iC=1:nCells  % for cells in zone
            % get points
            LOCATIONS{izone}(k+(1:nPoints),:)=...
                [ones(nPoints,1)*CRL(iC,:) relLoc ones(nPoints,1)*[zeros(1,3) t_released]];
            k=k+nPoints;
        end
    end

%    LOCATIONS= cell2list(LOCATIONS);

end

function relLoc = distribute(n,iFace)
    % relloc= distribute(n,iFace) --- compute relative locations of staring
    % points given n = [Nx,Ny,Nz] and iFace
    % TO 120528

    % accept n = instead of [Nx .. Nz]
    Nx=max(1,n(  1));
    Nz=max(1,n(end));
    
    % accept n = instead of [Nx Ny Nz]
    if length(n)>1,  Ny=n(2);
    else             Ny=n(1);
    end
    Ny=max(1,Ny);
    
    nFace = length(iFace);
    
    % distribute (relative coordinates)
    Ix = -1/(2*Nx) + (1:Nx)/Nx; 
    Iy = -1/(2*Ny) + (1:Ny)/Ny;
    Iz = -1/(2*Nz) + (1:Nz)/Nz;
    
    [IX,IY,IZ]=meshgrid(Ix,Iy,Iz);
    
    IX = IX(:)*ones(1,nFace);
    IY = IY(:)*ones(1,nFace);
    IZ = IZ(:)*ones(1,nFace);
    
    % allow iFace to be a vector for simultaneously defining points at more
    % than one iFace
    for i = 1:length(iFace)
        switch iFace(i)
            case 1, IX(:,i) = 0;
            case 2, IX(:,i) = 1;
            case 3, IY(:,i) = 0;
            case 4, IY(:,i) = 1;
            case 5, IZ(:,1) = 0;
            case 6, IZ(:,i) = 1;
            otherwise
                % nothing, i.e distribute across cell
        end
    end
    relLoc = unique([IX(:) IY(:) IZ(:)],'rows');
end
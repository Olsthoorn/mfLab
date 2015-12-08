function [xR yR zR] = particleRelCoords(iface,placement)
    % [xR yR zR] = particleRelCoords(iface,placement) -- compute the
    % relative coordinates [xR yR zR] for modpath particles within a cell
    % given iface and placement.
    % placement is the subdivision of the particles given an iface [0..6] and
    % the actual subsivsion, which is a 3 element vector if iface==0
    % or a 2 element vector in all other cases. Scalars may also be
    % used.
    % TO 130130
    
    placement(end:3) = placement(end);

    nL = placement(1);
    nR = placement(2);
    nC = placement(3);

    xRel=(-1+1/nC:2/nC:1-1/nC)/2;
    yRel=(-1+1/nR:2/nR:1-1/nR)/2;
    zRel=(-1+1/nL:2/nL:1-1/nL)/2;
    
    switch iface
        case 1
            xRel=-0.5;
        case 2
            xRel=+0.5;
        case 3
            yRel=-0.5;
        case 4
            yRel=+0.5;
        case 5
            zRel=-0.5;
        case 6
            zRel=+0.5;
    end
    
    % convert to full grid and shift so that relative coordinates vary from
    % zero to one withing the cell.
    [xR,yR,zR]=meshgrid(xRel,yRel,zRel);
end

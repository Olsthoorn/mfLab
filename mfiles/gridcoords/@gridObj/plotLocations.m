function [x,y,z,time] = relloc2global(o,locations,varargin)
    % gridObj/startLoc: h = gr.plotLocations(locations,varargin)
    % locations are points for or from MODPATH
    %
    % h is a handle to the plotted locations
    %
    % TO 120528
    
    %% j stands for column number
    jC = 1; jR=2; jL =3; jxRel=4; jyRel=5; jzRel=6; jtime=10; % jcCode=7; jrCode=8; jlCode=8;
    
    x = o.xGr(locations(:,jC)) + locations(:,jxRel).*o.dx(locations(:,jC));
    y = o.yGr(locations(:,jR)) + locations(:,jyRel).*o.dy(locations(:,jR));
    
    %% +1 at JL is because rel coordinate in z-direction counts from bottom of cell
    Idx = cellIndex(locations(:,jC), ...
                    locations(:,jR), ...
                    locations(:,jL + 1), gr.size);
    
    z = o.Z(Idx) + o.DZ(Idx) .* locations(:,jzRel);
    
    time = locations(:,10);
        
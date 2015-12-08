function modelLoc = relloc2model(o,locations)
% modelLoc = relloc2model(o,locations)
% locations are relative according to MODPATH:
%  [ Col Row Lay xRel yRel zRel 0 0 0 t_released]
%
% modelLoc = [x y z t_relased]
% locations and modelLoc may be a cell array with packages of locations
% i.e. with sublists or a list
% CellArray in = CellArra out
%
% h is a handle to the plotted locations
%
% TO 120528
    
    if iscell(locations)
        modelLoc{length(locations),1}={};
        
        for i=1:length(locations)
            modelLoc{i} = rel2mod(o,locations{i});
        end
    else
        modelLoc = rel2mod(o,locations);
    end
end

function modelLoc = rel2mod(o,locations)
% modelLocations = gr.rel2mod(relativeLocations)
% Convert relative locations (as list) to model locations (as list)
%
% TO 120528

%% j stands for column number
    jC = 1; jR=2; jL =3; jxRel=4; jyRel=5; jzRel=6; jtime=10; % jcCode=7; jrCode=8; jlCode=8;

    x = o.xGr(locations(:,jC)  )' + locations(:,jxRel).*o.dx(locations(:,jC))';
    y = o.yGr(locations(:,jR)+1)  + locations(:,jyRel).*o.dy(locations(:,jR));

    %% +1 at JL is because rel coordinate in z-direction counts from bottom of cell
    Idx = cellIndex(locations(:,jC), ...
                    locations(:,jR), ...
                    locations(:,jL), o.size);

    z = o.ZBlay(Idx) + o.DZlay(Idx) .* locations(:,jzRel);

    time = locations(:,jtime);
    
    modelLoc = [x y z time];

end
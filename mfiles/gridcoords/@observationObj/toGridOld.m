function o=toGrid(o,gr) 
%% obs = obs.toGrid(gr) places the observation point(s) in the grid
%
% TO 121118 121129

for i=length(o):-1:1
        
    [o(i).ix, o(i).iy, iz] = xyzindex([o(i).x(1),o(i).y(1),o(i).z(1)],gr);
    
    if any(isempty([o(i).ix,o(i).iy, iz]) | isnan([o(i).ix,o(i).iy,  iz]) | ~gr.isLay(iz))
        % o(i).remark='obs outside model'; % if you want to keep this obs
        warning('on','observationObj:toGrid:screenOutsideGrid');
        warning('observationObj:toGrid:screenOutsideGrid',...
            'removing obs(%d).[x,y,z]=[%g,%g,%g] because outside mesh or in a confining bed.\n',...
            o(i).nr,o(i).x(1),o(i).y(1),o(i).z(1));
         o(i)=[];
         warning('off','observationObj:toGrid:screenOutsideGrid');
         continue; % obs(iw)=[];
    else
        % Transfer the gridOlineObj contents to the obs screens
        o(i).x    = roundn(o(i).x,3);
        o(i).y    = roundn(o(i).y,3);
        o(i).z    = roundn(o(i).z,3);
        o(i).iLay = sum(gr.isLay(1:iz));
        o(i).idx  = cellIndex(o(i).ix,o(i).iy,o(i).iLay,gr.size);  
        o(i).LRC  = [o(i).iLay(:) o(i).iy(:) o(i).ix(:)];
    end
        
end

if isempty(o)
    error(['%s: all obs screens are at least partly outside the model !\n',...
          'Check z-coordinates of screens and then their x and y coordinates\n',...
          'Also verify the warnings ussued to the screen, telling which screen is out of the grid.']);
end

% Combine observations if they have the same id
   for i=numel(o)-1:-1
       % which one higher up has same number as this one ?
       % can at most be one because we started at top and work backwards
       j = find(o(i).id == [o(i+1:end).id]);
       
       if ~isempty(j)           
           j=j+i;           
            o(i).x    = [o(i).x    , o(j).x];
            o(i).y    = [o(i).y    , o(j).y];
            o(i).z    = [o(i).z    , o(j).z];
            o(i).ix   = [o(i).ix   , o(j).ix];
            o(i).iy   = [o(i).iy   , o(j).iy];
            o(i).iLay = [o(i).iLay , o(j).iLay];
            o(i).idx  = [o(i).idx  , o(j).idx];        
            o(i).LRC  = [o(i).LRC  ; o(j).LRC]; 
            
            o(j) = []; % remove o(i+1)
       end           
   end
end

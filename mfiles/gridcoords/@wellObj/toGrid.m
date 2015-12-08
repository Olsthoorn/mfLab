function o=toGrid(o,gr,HK) 
%% well = well.toGrid(gr,HK) places the MNW in the grid
% by computing indices and fraction of the discharge from the
% cells penetrated by the well screen.
%
% TO 121118

%% Deal with confining beds, that's tricky                            
LAYCBD = gr.LAYCBD;

%%
% if any LAYCBD, just use a new clean grid without LAYCBD
% expand HK to reflect the new grid and compute the well info as usual
if any(LAYCBD) % then we have confining beds
    grOld = gr;
    HK_old= HK;
    gr    = gridObj(grOld.xGr,grOld.yGr,grOld.Z); % new grid without confining beds
    HK    = gr.const(0);

    HK(:,:,grOld.isLay) = HK_old;
    
    isCbd = find(~grOld.isLay);
    HK(:,:,isCbd)= HK_old(:,:,isCbd-1);
    
    %% then continue as usual
end


for i=numel(o):-1:1

    if numel(o(i).x)<1 || numel(o(i).y)<1 || numel(o(i).z)<2
        error(['%s: wellObj(%d)\n',...
            ' has %d x coordinates, it must have at least 1;\n',...
            ' has %d y coordinates, it must have at least 1;\n',...
            ' has %d z coordinates, it must have at least 2 (top and bottom of well screen).\n'...
            'Check the wells in your call or the input table for the wells in your workbook.'],...
            mfilename,i,numel(o(i).x),numel(o(i).y),numel(o(i).z));
    end

    %see if the well has more than one coordinate for x and or y
    %then we have to cut this oblique line through the model and
    %set the paramters of the multi node well accordingly.
    % if it works well, it can also be used for ordinary wells
    
    if length(o(i).x)<numel(o(i).z)
        o(i).x(end+1:numel(o(i).z))=o(i).x(end); % extend x to match size of z
    end
    if length(o(i).y)<length(o(i).z)
        o(i).y(end+1:numel(o(i).z))=o(i).y(end); % extend y to match size of z
    end
    
    grdLines = gridLineObj(gr,[o(i).x(:) o(i).y(:) o(i).z(:)]);

    if isempty(grdLines)
        % o(i).remark='well outside model'; % if you want to keep this well
        warning('on','all');
        warning('wellObj:toGrid:screenOutsideGrid',...
            ['well(%d).[x,y,[z([1 end])]=[%g,%g, z =[%g,%g]] is (partly) outside the model mesh.\n',...
             'This well will be removed !'],o(i).nr,o(i).x(1),o(i).y(1),o(i).z(1),o(i).z(end));
         o(i)=[];
         continue; % well(iw)=[];
    else
        % Transfer the gridLineObj contents to the well screens
        o(i).x    = roundn([grdLines.x ],3);
        o(i).y    = roundn([grdLines.y ],3);
        o(i).z    = roundn([grdLines.z ],3);
        o(i).ix   = [grdLines.ix];
        o(i).iy   = [grdLines.iy];
        o(i).iLay = [grdLines.iz];
        o(i).idx  = [grdLines.idx];        
        o(i).ztop = max(gr.Zlay(o(i).iy(1),o(i).ix(1),:));
        o(i).LRC  = [o(i).iLay(:) o(i).iy(:) o(i).ix(:)];
        o(i).DZ   = gr.DZlay(o(i).idx);
        o(i).T    = o(i).DZ.*HK(o(i).idx);
        o(i).fQ   = o(i).T/sum(o(i).T)';
    end
        
end

if isempty(o)
    error(['%s: No wells remain in the model. This implies that all the\n',...
          'well screens are at least partly outside the model ! Therefore,\n',...
          'check z-coordinates of the screens and compare them with the grid layer elevations.\n',...
          'Then verify their x and y coordinates, and make sure that they are within the grid.\n',...
          'Also verify the warnings issued to the screen, telling which screen is out of the grid.'],...
          mfilename);
end

%% Translate the indices back to the original grid with the LAYCBD
% that is excluding the confining beds for the screens
if any(LAYCBD)
    % tranfer array, linking new with old layer indices
    transfer = cumsum([gr.isLay grOld.isLay]);
    
    for i=1:numel(o)
        % layers and coordinates in the newGrid that are valid in the old grid
        valid       = ismember(o(i).iLay,find(grOld.isLay));
        validCoords = [valid(:)'; valid(:)']; validCoords= validCoords(:)';
        
        % select valid coords
        o(i).x    = o(i).x(validCoords);
        o(i).y    = o(i).y(validCoords);
        o(i).z    = o(i).z(validCoords);
        
        % select valid indices and cell values
        o(i).ix   = o(i).ix(valid);
        o(i).iy   = o(i).iy(valid);
        o(i).iLay = transfer(o(i).iLay(valid),2)'; %right layer numbers
        o(i).idx  = o(i).idx(valid);   
        o(i).LRC  = o(i).LRC(valid,:);
        o(i).DZ   = o(i).DZ(valid);
        o(i).fQ   = o(i).fQ(valid);

        % renumber the layers
        o(i).LRC  = [o(i).iLay(:) o(i).iy(:) o(i).ix(:)];
        o(i).idx  = ((o(i).iLay(:)-1)*gr.Nxy + (o(i).ix(:)-1)*gr.Ny+o(i).iy(:))';
    end
end

% Combine screens if they have the same id (only MNW)
warning('on' ,'toGrid:wells_id_not_unique');
if numel([o.id]) ~= numel(unique([o.id]))
    warning('toGrid:wells_id_not_unique',...
        'wells with non unique id''s found, they will be combined, yielding fewer wells in the end.');
end
warning('off','toGrid:wells_id_not_unique');

warning('on','toGrid:removing_screen');
for i=numel(o)-1:-1:1       
   % which one higher up has same number as this one ?
   % can at most be one because we started at top and work backwards
   j = find(o(i).id == [o(i+1:end).id]);

      
   if ~isempty(j)
       warning('toGrid:removing_screen',...
           'well(%d) will be combined with well(%d), having the same id=%d',i,j+i,o(i).id);
       j=j+i;           
        o(i).x    = [o(i).x    , o(j).x];
        o(i).y    = [o(i).y    , o(j).y];
        o(i).z    = [o(i).z    , o(j).z];
        o(i).ix   = [o(i).ix   , o(j).ix];
        o(i).iy   = [o(i).iy   , o(j).iy];
        o(i).iLay = [o(i).iLay , o(j).iLay];
        o(i).idx  = [o(i).idx  , o(j).idx];        
        o(i).LRC  = [o(i).LRC  ; o(j).LRC]; 
        o(i).DZ   = [o(i).DZ   , o(j).DZ];
        o(i).fQ   = [o(i).fQ   , o(j).fQ];

        o(j) = []; % remove o(i+1)
   end           
end
warning('off','toGrid:removing_screen');



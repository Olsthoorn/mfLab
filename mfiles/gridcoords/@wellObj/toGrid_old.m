function o=toGrid(o,gr,HK) 
%% well = well.toGrid(gr,HK) places the well in the grid
% by computing indices and fraction of the discharge from the
% cells penetrated by the well screen.
%
% TO 120512

for i=length(o):-1:1

    [o(i).ix,o(i).iy] = xyzindex(o(i).x(:),o(i).y(:),gr.xGr,gr.yGr);

    if any(o(i).ix)>gr.Nx || any(o(i).iy)>gr.Ny || any(o(i).ix)<1 || any(o(i).iy)<1
        o(i)=[];
        continue;
    end
    
    if any(isnan(o(i).ix)) || any(isnan(o(i).iy))
        o(i).remark='well outside model';
        warning('mfsetwells:well:outsidemodel',...
            'well(%d).[x,y]=[%g,%g] is outside model mesh. This well will be removed !',o(i).nr,o(i).x,o(i).y);
         o(i)=[];
         continue; % well(iw)=[];
    else
        if size(gr.Zlay,2)==1, ix_=1; else ix_=o(i).ix; end
        if size(gr.Zlay,1)==1, iy_=1; else iy_=o(i).iy; end

        o(i).ztop=max(gr.Zlay(iy_,ix_,1));

        delta = 1e-3; % to make sure screen is within cell if it matches cell face exactly
        
        Lay1 = floor(interp1(XS(gr.Zlay(iy_(1),ix_(1),:)),1:(gr.Nlay+1),max(o(i).z)-delta)); % top    of well screen
        Lay2 = floor(interp1(XS(gr.Zlay(iy_(1),ix_(1),:)),1:(gr.Nlay+1),min(o(i).z)+delta)); % bottom of well screen
        Lay2 = min(Lay2,gr.Nlay);

        o(i).iLay=Lay1:Lay2; % non CBD layers penetrated by the screen

        if isnan(o(i).iLay)
            o(i).remark=[o(i).remark ', screen outside model at least parly'];
            warning('mfsetwells:screen:outsidemodel',...
                'well(%d).[x,y] [z1 z2]=[%g,%g] [%g %g] screen is (partly) outside model mesh. This part of entire well well will be removed !',...
                o(i).nr,o(i).x,o(i).y,o(i).z(1),o(i).z(end));
            o(i)=[];
            continue;
        end

        o(i).idx = cellIndex(o(i).ix(1),o(i).iy(1),o(i).iLay,size(HK));

        if ~isnan(o(i).idx)
            o(i).LRC = cellIndices(o(i).idx,size(HK),'LRC');
            o(i).DZ  = gr.DZlay(iy_,ix_,o(i).iLay);
            o(i).fQ  = (o(i).DZ(:).*HK(o(i).idx))/sum(o(i).DZ(:).*HK(o(i).idx))';
        end
    end
end

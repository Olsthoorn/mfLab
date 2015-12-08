function IFACE=IfaceIS(interface,Ny,Nx,NLAY,NSURF,TOP,BOTTOM)

IFACE(NLAY).ZETA=NaN(Ny,Nx,NSURF);
for i=1:NLAY
    for j=1:NSURF
        IFACE(i).ZETA(:,:,j)=interface(:,:,j);
        [rt,ct]=find(interface(:,:,j)>TOP(:,:,i));
        [rb,cb]=find(interface(:,:,j)<BOTTOM(:,:,i));
        for k=1:length(rt)
            IFACE(i).ZETA(rt(k),ct(k),j)=TOP(rt(k),ct(k),i);
        end
        for k=1:length(rb)
            IFACE(i).ZETA(rb(k),cb(k),j)=BOTTOM(rb(k),cb(k),i);
        end
    end
end


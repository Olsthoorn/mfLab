function [] = quiver3(o,B)


Qx = B(it).term{strmatchi('FLOWRIGHTFACE',B(it).label)};
Qy = B(it).term{strmatchi('FLOWFRONTFACE',B(it).label)};
Qz = B(it).term{strmatchi('FLOWLOWERFACE',B(it).label)};



Qx = cat(2,Qx(:,1,:),Qx); qx = (Qx(:,1:end-1,:)+Qx(:,2:end,:))/2;
Qy = cat(1,Qy(1,:,:),Qy); qy = (Qy(1:end-1,:,:)+Qy(2:end,:,:))/2;
Qz = cat(3,Qz(:,:,1),Qz); qz = (Qz(:,:,1:end-1)+Qz(:,:,2:end))/2;


quiver3(o.XM,o.YM,o.ZM,qx,qy,qz,1);
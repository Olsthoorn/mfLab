function showArray(gr,V,ix,iy,iz)
%SHOWARRAY shows an array with an outcut through cell ix,iy,iz
%
% Example:
%    showArray(gr,V,ix,iy,iz
%
% Used in:
%     mflab/examples/mf2k/Toth
%
% TO 120425

hold on
view(3)
    for iplane= 1:3

        [V1,V2]=getSurfaces(V,ix,iy,iz,iplane);

        switch iplane
            case 1
                surface(squeeze(gr.XMlay(end,:,:)), ...
                        squeeze(gr.YMlay(end,:,:)), ...
                        squeeze(gr.ZMlay(end,:,:)), ...
                        squeeze(V1(end,:,:)));
                surface(squeeze(gr.XMlay( iy,:,:)), ...
                        squeeze(gr.YMlay( iy,:,:)), ...
                        squeeze(gr.ZMlay( iy,:,:)), ...
                        squeeze(V2( iy,:,:)));
%             case 2
%                 surface(squeeze(gr.XMlay(:,end,:)), ...
%                         squeeze(gr.YMlay(:,end,:)), ...
%                         squeeze(gr.ZMlay(:,end,:)), ...
%                         squeeze(V1(:,end,:)));
%                 surface(squeeze(gr.XMlay(:, ix,:)), ...
%                         squeeze(gr.YMlay(:, ix,:)), ...
%                         squeeze(gr.ZMlay(:, ix,:)), ...
%                         squeeze(V2(:, ix,:)));
%             case 3
%                 surface(squeeze(gr.XMlay(:,:,  1)), ...
%                         squeeze(gr.YMlay(:,:,  1)), ...
%                         squeeze(gr.ZMlay(:,:,  1)), ...
%                         squeeze(V1(:,:, 1)));
%                 surface(squeeze(gr.XMlay(:,:, iz)), ...
%                         squeeze(gr.YMlay(:,:, iz)), ...
%                         squeeze(gr.ZMlay(:,:, iz)), ...
%                         squeeze(V2(:,:,iz)));
        end
    end
end
           
function [V1,V2]=getSurfaces(V,ix,iy,iz,iplane)

V1 = V;
V2 = NaN(size(V));

i=iplane;

switch i
   case 1, V1(  end , ix:end, 1:iz ) = NaN; V2(  iy , ix:end, 1:iz) = V( iy , ix:end, 1:iz);
   case 2, V1( 1:iy ,   end , 1:iz ) = NaN; V2( 1:iy,  ix , 1:iz )  = V(1:iy,   ix  , 1:iz);
   case 3, V1( 1:iy ,  1:ix ,   1  ) = NaN; V2( 1:iy, 1:ix,  iz  )  = V(1:iy,  1:ix ,  iz );
end
 
end
function IFace=getIface(Ny,Nx,NSURF,NLAY)
%GETIFACE obtain the interface elevation from the density values.
%
%
% first we will obtain them for the area where the data is provided (Zuid
% Holland) and afterwards we will put them into our model. Our area of
% interest is not fully contained in the given area but that will be solved
% later. In the meanwhile the elevation at the points which are not
% contained will be given a value NaN

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

pth='Z:\tolsthoorn On My Mac\GRWMODELS\NHI_and_AGV\mocdens3d';

load([pth filesep 'mocdense.mat']);

moc=rmfield(moc,'DZ');
moc=rmfield(moc,'SFLAG');
clear bas bcf

%% Remove all -1 values which are north sea water

moc.CHLOR(moc.CHLOR==-1)=max(moc.CHLOR(:));

%% Remove inversions by copying higher chlorides of overlying layer + DCHLOR
% this is primitive but soit

DCHLOR=10;


L1=zeros(size(moc.CHLOR));
L2=zeros(size(moc.CHLOR));
k=0;
while 1
    k=k+1;
    L1(:)=0; L2(:)=0;
    I=(moc.CHLOR(:,:,2:end)<=moc.CHLOR(:,:,1:end-1));
    fprintf('iteration %d, number of inversion points: %d\n',k,sum(I(:)));
    if any(I(:)),
        L1(:,:,1:end-1)=I;
        L2(:,:,2:end  )=I;
        moc.CHLOR(logical(L2))=moc.CHLOR(logical(L1))+DCHLOR;
    else
        break;
    end
end
clear L1 L2 I

%%
IFACEVALS=[300 1000 10000];

IFACE=zeros(moc.NROW,moc.NCOL,length(IFACEVALS));
CHLOR=NaN(moc.NROW,moc.NCOL,moc.NLAY+2);
CHLOR(:,:,2:end-1)=moc.CHLOR;
CHLOR(:,:,1  )=min(moc.CHLOR(:))   -1;  % actual minimum=5
CHLOR(:,:,end)=max(moc.CHLOR(:))+1000;  % actual maximum = 18700

Z=NaN(size(CHLOR)); Z(:,:,2:end-1)=moc.ZM; Z(:,:,1)=1000; Z(:,:,end)=-1000;
for ix=1:moc.NCOL
    fprintf('.')
    if rem(ix,50)==0, fprintf('%d/%d\n',ix,moc.NCOL); end
    for iy=1:moc.NROW
        IFACE(iy,ix,:)=interp1(squeeze(CHLOR(iy,ix,:)),squeeze(Z(iy,ix,:)),IFACEVALS);
    end
end

figure
for i=1:size(IFACE,3)
    subplot(3,1,i)
    contour(xm,ym,IFACE(:,:,i));
    title(sprintf('Interface %.0f mg/L',IFACEVALS(i)));
end

fclose('all');

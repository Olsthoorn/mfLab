function HFB = setHFB(xm,ym,ILay,xv,yv,c)
%SETHFB geneate the array required for HFPB processes.
%
% USAGE:
%    HFB = setHFV(Nx,Ny,ILay,xv,yv,c)
%
% Generate the array for necessary for the horiontal fow barrier process
% for layers ILay in the grid specified by xm,ym, horizontal flow boundary
% contour specified by the vertices xv yv and the resistance across this
% boundary given by c [time].
% ILay may contain the indices of more than one layer
% c may be of equal length as ILay, in which case c specifies the
% resistances for each layer separately. c will be expanded as necessary in
% case its length is smaller than that of ILay
% The resistance equals the reciprocal of the Hydrch variable in the
% modflow manuals. c = 1/Hydchr = d/k where d is the thickness of the
% barrier and k its horizontal conductivity.
%
% TO 101020
%
% Copyright 2009-2010 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later


%% in case xm is a row vector and ym a column vector make full 2D matrices   
nx=size(xm,2);
ny=size(ym,1);

if size(xm,1)==1, xm=ones(ny,1)*xm; end
if size(ym,2)==1; ym=ym*ones(1,nx); end

%% proceed

% figure
% plotgrid(xm(1,:),ym(:,1));

IN=inpolygon(xm,ym,xv,yv);
% plot(xm(IN),ym(IN),'k.');


EW=[diff(IN,1,2)  zeros(ny,1)];
NS=[diff(IN,1,1); zeros(1,nx)];

LRC_EW=cellIndices(find(EW~=0),[ny,nx],'LRC');
LRC_NS=cellIndices(find(NS~=0),[ny,nx],'LRC');

% for i=1:size(LRC_EW,1)
%     plot(xm(LRC_EW(i,2),LRC_EW(i,3)  ),ym(LRC_EW(i,2),LRC_EW(i,3)  ),'g.');
%     plot(xm(LRC_EW(i,2),LRC_EW(i,3)+1),ym(LRC_EW(i,2),LRC_EW(i,3)+1),'m.');
% end
% 
% for i=1:size(LRC_NS,1);
%     plot(xm(LRC_NS(i,2)  ,LRC_NS(i,3)),ym(LRC_NS(i,2)  ,LRC_NS(i,3)),'b.');
%     plot(xm(LRC_NS(i,2)+1,LRC_NS(i,3)),ym(LRC_NS(i,2)+1,LRC_NS(i,3)),'r.');
% end

% plotgrid(xm(1,:),ym(:,1));
% plot(xm(IEW),ym(IEW),'r.')
% plot(xm(INS),ym(INS),'b.');

hfb = [LRC_EW(:,2) LRC_EW(:,3) LRC_EW(:,2)    LRC_EW(:,3)+1; ...
       LRC_NS(:,2) LRC_NS(:,3) LRC_NS(:,2)+1, LRC_NS(:,3)    ...
       ];
   
u=ones(size(hfb(:,1)));

if length(c)<length(ILay)
    c(end:length(ILay))=c(end);
end 

HFB=[];
for i=1:length(ILay)
    HFB=[HFB; ...
        ILay(i)*u hfb c(ILay(i))*u ...
        ];
end

 

function Z=mf_Z_extend(Z,xGr,yGr)
%MF_Z_EXTEND extends x and y of Z to nodes instead of cell centers
%
% Example:
%    ZGR=mf_Z_extend(Z[,xGr[,yGr]])
%
% Only used in mflab/examples/mf2007/Khettarra_Jorf/mf_adaptUnderConstruction
%
% ToDo: This may easily be replaced by Matlab's bsxfun()
%
% See also: modelsize3 gridObj bsxfun
%
% TO 100602


if nargin<3, yGr=ones(size(Z(:,1,1))); else xGr=xGr(:)'; end
if nargin<2, xGr=ones(size(Z(1,:,1))); else yGr=yGr(:);  end

if size(Z,1)<size(yGr,1)
    ZGR=NaN(size(yGr,1),size(Z,2),size(Z,3));
    ZGR([1 end],:,:)=Z([1 end],:,:);
    ZGR(2:end-1,:,:)=0.5*(Z(1:end-1,:,:)+Z(2:end,:,:));
    Z=ZGR;
end
if size(Z,2)<size(xGr,2)
    ZGR=NaN(size(Z,1),size(xGr,2),size(Z,3));
    ZGR(:,[1 end],:)=Z(:,[1 end],:);
    ZGR(:,2:end-1,:)=0.5*(Z(:,1:end-1,:)+Z(:,2:end,:));
    Z=ZGR;
end

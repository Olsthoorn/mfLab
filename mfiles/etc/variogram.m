function [s,h]=variogram(xm,ym,X)
%VARIOGRAM compute a variogram from the data
%
% No idea what this exactly is
%
% TO 090314

% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later


N=length(X(:)); o=ones(N,1);

[x,y]=meshgrid(xm(:),ym(:)');

R=sqrt((y(:)*o'-o*y(:)').^2+(x(:)*o'-o*x(:)').^2);
RM=mean(R(:));
dX=X(:)*o'-o*X(:)';

bins=0:RM/20:RM; Nbins=length(bins)-1;

h=0.5*(bins(1:end-1)+bins(2:end));
s=zeros(1,Nbins);
for ib=1:Nbins
    s(ib)=var(dX(R>bins(ib) & R<=bins(ib+1)));
end



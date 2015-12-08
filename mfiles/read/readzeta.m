function zeta = readzeta(filename,nrow,ncol,nlay,npln,nstp,nprn)
%READZETA reads interfaces for Salt Water Intrusion package (SWI)
%
% Example:
%    zeta = Readzeta(filename,nrow,ncol,nlay,npln,nstp,nprn)
%    returns zeta(npln+2,nrow,ncol,nlay,nstp/nprn)
%
% TO 100101

fid = fopen(filename);
out = fread(fid,nrow*ncol*nlay*(npln+2)*nstp/nprn,'float');
out = reshape(out,ncol,nrow,nlay,npln+2,nstp/nprn);
zeta = permute(out,[4 2 1 3 5]);

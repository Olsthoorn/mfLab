function [BCN names]= getNHIBCN(codes,basename,sheetNm,columnHdr,filesdir,Ix,Iy)
%GETNHIBCN gets boundary conditions/stresses from files RIV GHB DRN
%
% Example:
%    BCN = getNHIBCN(type,codes,basename,sheetNm,columnNr)
%  
% TO 120430

[fname names         ] = getNHIfileNm(basename,sheetNm,columnHdr,codes{1});
STAGE = getNHIASC(fullfile(filesdir,fname),Ix,Iy);

[fname names(end+1,:)]= getNHIfileNm(basename,sheetNm,columnHdr,codes{2}); 
C     = getNHIASC(fullfile(filesdir,fname),Ix,Iy);
%% Note STAGE and C have the size of the subgrid (1:length(Iy),1:length(Iy))

I     = find(~isnan(STAGE));
LRC   = cellIndices(I,size(C),'LRC');
u     = ones(numel(I),1);

if length(codes)<3
    BCN    = [u LRC ...
                reshape(STAGE(I),[numel(I),1]), ...
                reshape(C(I),    [numel(I),1])];
else
    fname  = getNHIfileNm(basename,sheetNm,columnHdr,codes{3}); 
    BOTTOM = getNHIASC(fullfile(filesdir,fname),Ix,Iy);
    BCN    = [u LRC ...
                reshape(STAGE(I) ,[numel(I),1]),...
                reshape(C(I)     ,[numel(I),1]),...
                reshape(BOTTOM(I),[numel(I),1])];
end

%% BCN has the L R C of the subgrid (the model grid)

function BCN =STRTHD2BCN(BCN,type,HDS)
%% BCN = connect2STRTHD(BCN,'type',H)
% BCN is RIV, DRN, GHB ...
% H is the head struct as it is read in with readDat([basename '.HDS']);
% The head in the BCN will be replaced with the SRTHD for the same stress
% period
% TO 120614

depth = 1.5; 

k=0;
for it=1:length(HDS)
    I=find(BCN(k+1:end,1)==it);
    Idx = cellIndex(BCN(k+I,4),BCN(k+I,3),BCN(k+I,2),size(HDS(1).values));
    switch type
        case {'DRN','GHB'}
            BCN(k+I,5)=HDS(it).values(Idx);
        case 'RIV'
            BCN(k+I,5)=HDS(it).values(Idx);
            BCN(k+I,7)=HDS(it).values(Idx)-depth;
        case 'CHD',
            BCN(k+I,5)=HDS(it).values(Idx);
            BCN(k+I,5)=HDS(it).values(Idx);
        otherwise
            error('%s: unknown boundary type %s',mfilename,type);
    end
end

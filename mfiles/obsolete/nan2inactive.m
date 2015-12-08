%function nan2inactive
%% remove nans for cells that are inactive and set the corresponding IBOUND cells equal to 0

% Just get all arrays that have size of IBOUND

IBOUND(isnan(STRTHD))=0;

STRTHD(isnan(STRTHD))= 0;
TRAN  (isnan(TRAN )) = 0;
VCONT (isnan(VCONT)) = 0;
RECH  (isnan(RECH )) = 0;


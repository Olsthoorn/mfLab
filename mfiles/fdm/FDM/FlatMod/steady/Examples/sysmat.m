function A=sysmat(kD,c)
% functie sysmat(kD,c) stelt de systeem-matrix samen.
% kD en c zijn kolomvectoren
n=length(kD);
a=1./(kD.*c);
b=1./(kD.*[c(2:n);inf]);
A=diag(a+b)-diag(a(2:n),-1)-diag(b(1:n-1),1);

function [Z,dy]=nicerange(Y,N)
% gets a nice equidistance range over y in N equal steps
% TO 000528

r=[1;2;5]*10.^[-5:5];
M=max(Y(:));
m=min(Y(:));
dy=(M-m)/max(2,fix(N));
dy=r(max(find(r(:)<=dy)));
n1=fix(min(Y(:))/dy);
n2=fix(max(Y(:))/dy);
Z=dy*[n1:n2];
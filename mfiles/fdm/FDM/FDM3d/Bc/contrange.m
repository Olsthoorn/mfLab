function [yrange,dy]=contrange(Y,n);
%[yrange,dy]=contrange(Y,n);
% create a decent range yrange and contour stepsize dy for contouring, using Y and n as number of rangesteps
% TO 000120

d=[0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000,2000,5000];

m=floor(min(Y(:)));
M= ceil(max(Y(:)));
dy=(M-m)/n;
dy=d(max(find(dy>=d)));
yrange=[m:dy:M];

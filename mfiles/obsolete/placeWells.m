function [xp yp]=plaatsputten(xy,N)
%function [xp yp]=plaatsput(xy,N)
% plaats N putten langs geplotte lijn [xp yp]
% TO 091123

r=[0:1/(N-1):1];
xp=xy(1,1)+r*(xy(end,1)-xy(1,1)); xp=xp(:);
yp=interp1(xy(:,1),xy(:,2),xp);   yp=yp(:);

fprintf('%12g\t%12g\n',[xp(:) yp(:)]');


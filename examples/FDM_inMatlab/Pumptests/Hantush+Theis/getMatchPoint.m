function getMatchPoint()
%GETMATCHPOINT  analyze classic pumping test interactively on screen
%
% USAGE getMatchPoint() % after moving points on screen.
%       function is invoked by amimator_loglog at button up,
%       i.e. button release of the data points moved across the screen.
%
% Example:
%    getMatchPoint() 
%
% what is does
%    1) gets the new and old matchpoints
%    2) computes the transmissvity and storage coefficient
%      assuming s/Q on vertical log axis and t/R^2 on
%      horizontal log axis.
%
% See also: Pumptest animator_loglog
%
% TO 120116

 h0=findobj('Tag',        'MatchPoint');
 h1=findobj('DisplayName','MatchPoint');

 x0=get(h0,'xData');
 y0=get(h0,'yData');
 
 x1=get(h1,'xData');
 y1=get(h1,'yData');
 
 F1=y1/y0;
 F2=x1/x0;

 T=F1/(4*pi);
 S=4*T/F2;
 
 fprintf('Transmissivity = %.4g\n'  ,T);
 fprintf('Storage coef   = %.4g\n\n',S);
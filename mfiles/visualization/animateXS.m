%ANIMATEXS script to get an animateObj and to run an animation with movie creation
%
% Example:
%    animate = animateObj(basename,gr,{'Salinity','heads','budget'});
%    animate.concXS(well);
%
%
% See also: animateObj
%
% TO 120420 130407

load name;
load(basename);

animate = animateObj(basename,gr,{'Salinity','heads','budget'});

animate.concXS(well);

        
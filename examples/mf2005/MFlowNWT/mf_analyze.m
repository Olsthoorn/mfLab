% MF_ANALYZE: Visualize and interpret model output
%
% USAGE:
%   mf_analyze
%
% MFlow NWT - Newton a Newton Formulation for MODFLOW 2005
% Richard G. Niswonger, Sorab Panday, USGS 2011, Chapter 37 of section A
% groundwater Book 6: Modeling Techniques.
%
% Analyzing output of the model
%
% TO 130425

clear variables; close all;

load name;
load(basename);

animate = animateObj(basename,'heads');

IRow=1;

animate.headXS(gr,IRow)

%% plots on cross section

set(gca,'nextplot','add','ylim',gr.zGr([end 1]));
xlabel('r [m]'); ylabel('elevation');
title('head in all layers for all stress periods');
for it=1:numel(animate.H)
    plot(gr.xm,squeeze(animate.H(it).values(1,:,:)),mf_color(it));
end

% Plot the data published in 
data = [  62.5 62.58
        1062.5 43.59
        2062.5 37.07
        3062.5 32.47];
plot(data(:,1),data(:,2),'ro','markerfacecolor','r');

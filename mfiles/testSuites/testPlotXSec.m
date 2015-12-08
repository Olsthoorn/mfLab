%TESTPLOTXSEC tests plotXSec and plotYSec for many input variants
%
% Example:
%    h=gr.plotXSec(varargin) -- plots XSections along xaxis,
%
%   through all Jrow of the mesh. These cross secitons are plotted on their
%   own figure/axis or subplot axis.
%
% Aquitards are filled with a sequence of colors. This is automatically the
% case of confining beds are present, that is, any gr,LAYCBD ~=0.
%
% Aquifers may also be filled if ~isempty(ILay).
%
% 
% USAGE
%  [hlay hcbd] = gr.plotXSec(Jrow,Ilay,options)
%  [hlay hcbd] = gr.plotXSec(Jrow,Ilay,Icbd,options)
%  [hlay hcbd] = gr.plotXSec(Jrow,'all',plotoptions)
%  [hlay hcbd] = gr.plotXSec(Jrow,'lay',Ilay,options)
%  [hlay hcbd] = gr.plotXSec(Jrow,'lay',LAYVAR,'cbd',CBDVAR,'clim',clim,options)
%
% Jrow: (first numeric input)  vector of the cross sections and or layers to be plotted.
% Ilay: (second numeric input) indices or logical vector telling which model layers to plot
% Icbd: (third numeric input)  indices or logical vector telling which confining beds to plot 
%
% options:
%   'figure'  --> new figure requested
%   'all'     --> color all layers and confining beds
%   'lines' ,cllorSpec  --> plot both hlines and vlines (cel boundaries)
%   'hlines',colorSpec  --> plot lines between layer in specified color
%   'vlines',colorSpec  --> plots vertical lines between cells in specified color
%          NOTE that colorSpec may be also 'on' and 'off' or 'none'
%   'lay',3Darray  --> to use for coloring XSec
%   'cbd',3Darray  --> to use for coloring XSec
%   'smooth' --> plot layer interfaces as contniuous lines
%   'stair'  --> plot layer interfaces as stair lines (showing actual cell bottoms and tops)
% plotOptions
%
% plotOptions (any valid options for Matlab's fill and patch functions)
%   'edgecolor',colorSpec
%   'faceColor',colorSpec
%   'edgeAlpha',fraction
%   'faceAlpha',fraction
%   'lineWidth',value
%   'color',colorSpec
%   plus any options that are accepted by matlab's fill and patch
%
% NOTE
% Ilay and Icbd may be vectors of layer indices or logical vectors, i.e. the
% result of a logical expression. The latter allows parameterizing the layers
% to be filled. For instance, by demanding that its conductivity is less than
% some prescribed value.
% For instance a criterion like this
%    max(max(HK(:,:,:),[],1),[],2)<0.03
% is such a logical vector selecting only the layers whose max HK<0.03.
%
%
% SEE ALSO: gr.plotYSec | gr.plotGrid | gr.mesh | gr.plotXS | gr.plotYS
%
% TO 120501 120531 130208

clear; close all;


%% Generate a random model to show cross sections.

xGr = -1000:50:1000; Nx = numel(xGr)-1;
yGr = -1000:50:1000; Ny = numel(yGr)-1;
zGr = -100:10:0;     Nz = numel(zGr)-1;
gr = gridObj(xGr,yGr,zGr);

k = rand(Ny,Nx,Nz);
FH = NaN(  size(k)); FH(:,:,end) = 0;
FQ = zeros(size(k)); FQ(:,:,  1) = 1000;
LAYCBD=[0 1 0 0 0 6 0];
%LAYCBD=[1 3 6 7];

Z = fdm3(gr.xGr,gr.yGr,gr.Z,k,k,k,FH,FQ)-100;
Z = Z - repmat(Z(:,:,round(gr.Nz/2)),[1,1,gr.Nz]);

gr = gridObj(xGr,yGr,Z,LAYCBD);

% generate a variable conducivity field
HK   = gr.ZMlay; for iLay=1:gr.Nlay, HK(:,:,iLay) = abs(HK(:,:,iLay)-mean(mean(HK(:,:,iLay)))); end
VKCB = HK(:,:,1:gr.Ncbd)/1000;

%% EXAMPLEs / tests:

TEST = 14;

fprintf('\n\nTesting plotXSec:\n\n');

for test = TEST
    
    close all

    switch test
        case 1
            gr.plotXSec([1 3 5 500 -3 2]);
        case 2
            gr.plotXSec(5,'fig');
        case 3
            gr.plotXSec(5,3,'fig');
            gr.plotXSec(5,[],[1 3],'fig');
            gr.plotXSec(5,[2 5],[],'fig');
            gr.plotXSec(5,[],[],'fig');
        case 4
            gr.plotXSec(5,[],[]);
        case 5
            gr.plotXSec(5,'lay',HK,'cbd',VKCB);
        case 6
            gr.plotXSec(1,[5 7 9],'lay',HK,'cbd',VKCB,'hlines','g','vlines','b','linewidth',2,'facecolor','g');
            gr.plotXSec('fig',1,'smooth',[5 7 9],'lay',HK,'cbd',VKCB,'hlines','g','vlines','b','linewidth',2,'facecolor','g');
        case 7
            gr.plotXSec(1,[5 7 9],'lay',HK,'cbd',VKCB,'hlines','m','vlines','b','linewidth',2,'facecolor','g');
        case 8
            gr.plotXSec(1,[5 7 9],'smooth','hlines','y','vlines','b','linewidth',2,'facecolor','g');
        case 9
            gr.plotXSec(1,[5 7 9],'lay',log10(HK),'hlines','k','faceAlpha',0.25);
        case 10
            gr.plotXSec(1:10:gr.Ny,'all','figure','vlines','g','smooth');
        case 11
            gr.plotXSec(1,mean(mean(HK,2),1)<15,'stair')
        case 12
            gr.plotXSec(1:5:gr.yGr,'hlines',[0.8 0.8 0.8],'title','test plotXSec','smooth');
            gr.plotXSec(1:5:gr.yGr,'hlines',[0.8 0.8 0.8],'title','test plotXSec');
    end

end


%% testing plotYSec
fprintf('\n\nTesting plotYSec:\n\n');

for test = TEST
    
    close all

    switch test
        case 1
            gr.plotYSec([1 3 5 500 -3 2]);
        case 2
            gr.plotYSec(5,'fig');
        case 3
            gr.plotYSec(5,3,'fig');
            gr.plotYSec(5,[],[1 3],'fig');
            gr.plotYSec(5,[2 5],[],'fig');
            gr.plotYSec(5,[],[],'fig');
        case 4
            gr.plotYSec(5,[],[]);
        case 5
            gr.plotYSec(5,'lay',HK,'cbd',VKCB);
        case 6
            gr.plotYSec(1,[5 7 9],'lay',HK,'cbd',VKCB,'hlines','g','vlines','b','linewidth',2,'facecolor','g');
            gr.plotYSec('fig',1,'smooth',[5 7 9],'lay',HK,'cbd',VKCB,'hlines','g','vlines','b','linewidth',2,'facecolor','g');
        case 7
            gr.plotYSec(1,[5 7 9],'lay',HK,'cbd',VKCB,'hlines','m','vlines','b','linewidth',2,'facecolor','g');
        case 8
            gr.plotYSec(1,[5 7 9],'smooth','hlines','y','vlines','b','linewidth',2,'facecolor','g');
        case 9
            gr.plotYSec(1,[5 7 9],'lay',log10(HK),'hlines','k','faceAlpha',0.25);
        case 10
            gr.plotYSec(1:10:gr.Ny,'all','figure','vlines','g','smooth');
        case 11
            gr.plotYSec(1,mean(mean(HK,2),1)<15,'stair')
        case 12
            gr.plotYSec(1:5:gr.yGr,'hlines',[0.8 0.8 0.8],'title','test plotYSec','smooth');
            gr.plotYSec(1:5:gr.yGr,'hlines',[0.8 0.8 0.8],'title','test plotYSec');
        case 13
            gr.plotYSec(1,'fig','lines','on','title','lines on plotYSec','smooth');
            gr.plotXSec(2,'fig','lines','none','title','lines none plotXSec','smooth');
            gr.plotYSec(3,'fig','lines','off','title','lines off plotYSec','smooth');
            
            gr.plotYSec(1,'fig','hlines','on','title','hlines on plotYSec','smooth');
            gr.plotYSec(2,'fig','hlines','none','title','hlines none plotYSec','smooth');
            gr.plotYSec(3,'fig','hlines','off','title','hlines off plotYSec','smooth');
            
            gr.plotYSec(1,'fig','vlines','on','title','vlines on plotYSec','smooth');
            gr.plotYSec(2,'fig','vlines','off','title','vlines ooff plotYSec','smooth');
            gr.plotYSec(3,'fig','vlines','none','title','vlines none plotYSec','smooth');
            gr.plotYSec(4,'fig','hlines','none','vlines','off','title','hlines none vlines off test plotYSec','smooth');
            gr.plotYSec(5,'fig','hlines',[0.8 0.8 0.8],'vlines','on','title','hlines grey vlines on test plotYSec');            
        case 14
            gr.plotYSec(1,'fig','lines','on','title','lines on plotYSec','stair');
            gr.plotXSec(2,'fig','lines','none','title','lines none plotXSec','stair');
            gr.plotYSec(3,'fig','lines','off','title','lines off plotYSec','stair');
            
            gr.plotYSec(1,'fig','hlines','on','title','hlines on plotYSec','stair');
            gr.plotYSec(2,'fig','hlines','none','title','hlines none plotYSec','stair');
            gr.plotYSec(3,'fig','hlines','off','title','hlines off plotYSec','stair');
            
            gr.plotYSec(1,'fig',[],'vlines','on','title','vlines on plotYSec','stair');
            gr.plotYSec(2,'fig',[],'vlines','off','title','vlines ooff plotYSec','stair');
            gr.plotYSec(3,'fig',[],'vlines','none','title','vlines none plotYSec','stair');
            gr.plotYSec(4,'fig',[],'hlines','none','vlines','off','title','hlines none vlines off test plotYSec','stair');
            gr.plotYSec(5,'fig',[],'hlines',[0.8 0.8 0.8],'vlines','on','stair','title','hlines grey vlines on test plotYSec');            
    end

end

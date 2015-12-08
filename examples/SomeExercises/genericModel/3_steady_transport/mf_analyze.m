%% Visualization of steady model results with particle tracking

%% Closing existing graphs and variables
close all;
clear variables;

%% The options and steps will be outlined here
% Notice that there are many more visualization options. One of them are
% animations and video's that are generally used in combation with
% transient flow and with mass transport simulations. Here we only have a
% steady-steate model.

%% First step: retrieve the basename of this model
load name       % retrieve basename stored in file name.mat
load(basename)  % get model arrays that were saved by mf_setup

%% Get the heads and the budget info.
H = readDat([basename '.HDS']);  % read computed heads

Tracer{1} = readMT3D('MT3D001.UCN');
Tracer{2} = readMT3D('MT3D002.UCN');

crange{1} = ContourRange(Tracer{1},50);
crange{2} = ContourRange(Tracer{2},50);

B = readBud([basename '.BGT']);  % read ceomputed budgets

%% Ways to show the data

%% Show zoneBudget
%
% generate a 3D zone array with zones 1, 2 and 3 corresponding
% to layers 1, 2, and 3. Then request the zoneBudget
zoneArray = gr.const([1 2 3]);   % generate the zoneArray
zonebudget(B(end),zoneArray,1);       % request zonebudget for layer 1
zonebudget(B(end),zoneArray,2);       % same for layer 2
zonebudget(B(end),zoneArray,3);       % same for layer 3
zonebudget(B(end),zoneArray,[1 2 3]); % same for all layers combined

%% Plotting conc in cross sections
for iComp = 1:2
    % along x-axis
    gr.plotXSec(1,'figure','xsec','title','Conc Tracer%d along y-axis','fontsize',14,'all','lay',Tracer{iComp}(end).values);
    hb = colorbar; set(get(hb,'title'),'string','Conc');

    % along y-axis
    gr.plotYSec(1,'figure','ysec','title','Conc Tracer%d along x-axis','fontsize',14,'all','lay',Tracer{iComp}(end).values);
    hb = colorbar;  set(get(hb,'title'),'string','Conc');
end

%% Plotting vertically averated conc contours of all aquifers
Ia{1} = gr.zmlay>-150 & gr.zmlay< Inf;
Ia{2} = gr.zmlay>-200 & gr.zmlay<-150;
Ia{3} = gr.zmlay>-450 & gr.zmlay<-350;

NCOMP = length(STCONC);

figure('name','contours tracer concentration','position',screenPos(0.6)); % fig = 60% of screen
for iComp=NCOMP:-1:1
    for ia=numel(Ia):-1:1
        j = (ia-1)*NCOMP+iComp;
        ax(ia,iComp) = subplot(numel(Ia),NCOMP,j,'nextplot','add');  % axis to plot layer
        set(ax(ia,iComp),'clim',[0 crange{iComp}(end)]);
    end
end

time = [Tracer{1}.time];

ht = NaN(numel(Ia),NCOMP);
hc = cell(numel(Ia),NCOMP);

st1 = 'average conc of Tracer%%d, aquifer%%%%d, time = %.4g d';

%%
viewDir = 'z';

for it=1:length(time)
    
    st2 = sprintf(st1,time(it));
    
    for iComp = 1:NCOMP

        st3 = sprintf(st2,iComp);
        
        hrange = ContourRange(H,50); % get a useful set of contour elevations

        % Plot the contours for the three layers
        for ia=1:numel(Ia)
            
            st4 = sprintf(st3,ia);
            
            if it==1
                xlabel(ax(ia,iComp),'x [m]');
                ylabel(ax(ia,iComp),'y [m]');  % axis labels
                ht(ia,iComp) = title(ax(ia,iComp),st4);      % title
                h = colorbar('peer',ax(ia,iComp)); set(get(h,'title'),'string','Conc');    % colorbar for heads

                % contours of aquifer-average conc
                switch viewDir
                    case 'x'
                        [~,hc{ia,iComp}] = contourf(ax(ia,iComp),gr.ym,gr.zm(:),YS(mean(Tracer{iComp}(it).values,2)),crange{iComp});
                    case 'y'
                        [~,hc{ia,iComp}] = contourf(ax(ia,iComp),gr.xm,gr.zm(:),XS(mean(Tracer{iComp}(it).values,1)),crange{iComp});
                    case 'z'
                        [~,hc{ia,iComp}] = contourf(ax(ia,iComp),gr.xm,gr.ym,mean(Tracer{iComp}(it).values(:,:,Ia{ia}),3),crange{iComp});
                    otherwise
                        error('%s: viewDir must be one of ''x'' ''y'' ''z'', not %s',mfilename,viewDir);
                end
                % plot well locations of wells pertaining to each of the layers
                % plotting is done using method plotXY of wellObj
                % selection of wells is done using [well.iLay]==ia
                lw = false(size(well));
                for iw=numel(well):-1:1
                    lw(iw)=any(ismember(well(iw).iLay,Ia{ia}));
                end
                
                well(lw).plotXY(ax(ia,iComp),'marker','o','markerEdgeColor','r');
            else
                set(ht(ia,iComp),'string',st4);
                switch viewDir
                    case 'x'
                        set(hc{ia,iComp},'zData',YS(mean(Tracer{iComp}(it).values,2)));                
                    case 'y'
                        set(hc{ia,iComp},'zData',XS(mean(Tracer{iComp}(it).values,1)));                
                    case 'z'
                        set(hc{ia,iComp},'zData',mean(Tracer{iComp}(it).values(:,:,Ia{ia}),3));                
                    otherwise
                        % viewDir already checked
                end
            end
        end
        drawnow();
    end
end

%% TO 130614
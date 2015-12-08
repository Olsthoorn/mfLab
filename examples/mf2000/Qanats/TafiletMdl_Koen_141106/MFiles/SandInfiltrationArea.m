%% Sand infiltration area
% used to get the sand infiltration areas

if 0
    %% This part has been done, and is saved to excel!!
    figure; hold on;
    axProps = {'nextplot','add','xGrid','on','yGrid','on'}; % default figure porporties
    XLIM = gr.xGr([1 end]);
    YLIM = gr.yGr([end 1]);
    axes(axProps{:},'xlim',XLIM,'ylim',YLIM);  % create axis, using props and limits
    xlabel('x UTM [m]'); ylabel('y UTM [m]');
    title('Heads in first layer at end of simulation period');
    A = imread('modelAreaImage.png'); image(gr.xGr,gr.yGr,A);
    title('Choose Oukhit Area'); axis equal;
    set(gcf, 'renderer', 'zbuffer');
    modelArea.plot
    riversQ.plot
    khettaras.plot

    [xx,yy] = ginput;
    [x,y] = utm2wgs(xx,yy,'30N');
    ExcelInfo = [x(:) y(:) xx(:) yy(:)];
end


%% Visualizing the output of EX2 of the Modflow 2000 manual.
% TO 091011 091129 120413 140112
 
%% Get data
load('name.mat') % get basename stored in file name.mat
load(basename);  % having retrieved baasename load the data in basename.mat
load underneath  % to get gr object

defaults = {'nextplot','add','fontsize',14,'xGrid','on','yGrid','on','fontSize',10};

% Read head file
H=readDat([basename  '.hds']);
B=readBud([basename, '.BGT']);

if gr.Ny==1
    B=mf_Psi(B);

    %% Plot isolines of head
    figure('pos',screenPos(0.6));

    % Create axes at their desired position
    ax(1)=axes(defaults{:},'ylim',[gr.zGr(end) max(H(end).values(:))]);
    xlabel('x [m]'); ylabel('elevation [m]');
    title('heads and stream lines, test fixed-gradient boundary conditions');

    plot(ax(1),gr.xm,squeeze(XS(H(end).values)));    
    contour(ax(1),gr.xp,gr.zp,B(end).Psi,50);
else
    %% Plot isolines of head
    figure('pos',screenPos(0.6));

    % Create axes at their desired position
    ax(1)=axes(defaults{:},'ylim',[gr.zGr(end) max(H(end).values(:))]);
    xlabel('x [m]'); ylabel('elevation [m]');
    title('heads and stream lines, test fixed-gradient boundary conditions');
    
    for iL=1:gr.Nz
        for iR=1:gr.Ny
            if ~any(isnan(H(end).values(iR,:,iL)))
                plot(ax(1),gr.xm,squeeze(XS(H(end).values(iR,:,iL))),[mf_color(iR),mf_linetype(iL),mf_marker(iR)]);
            end
        end
    end
end
%% gradient
grad = diff(H(end).values(:,:,:),1,2)./diff(gr.XM,1,2);
grad = grad(:,[1 end],:);
fprintf('\n\nwas set in mf_adapt : %10s %10.3g%10.3g\n','gradL gradR',gradientN,gradientS);
fprintf('obtained from model:\n %10s %10s\n','Left','Right');
for iR=1:gr.Ny
    fprintf('iRow = %d\n',iR);
    for iz=1:size(grad,3)
        fprintf(' %10.2g',grad(1,[1 end],iz));
        fprintf('\n');
    end
end
%% Use zone budget to get budget overview

zonebudget(B);

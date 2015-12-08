function showBox(o,tts)
%% Model.verify(tts) -- Verify model: check various aspects of the AMWADU model for plausibility
% tts = optional string for the figure header

gr = o.grid();

if ~all(gr.LAYCBD==0),
    warning('ModelObj:verify:cbdNotRemoved',...
        '%s: only works when Ncbd==0 (no confining beds). Use plotMesh or first Model.removeCBD',mfilename);
    return;
end

Iy = (1:gr.Ny)';

%%

figure('name','verfify','position',screenPos(0.75));
ax = axes('nextplot','add');
xlabel(ax,'x [m]'); ylabel(ax,'y [m]'); zlabel(ax,'z [m]');
view(3);

%% Title
if ~exist('tts','var'),
    tts='Verifying the Model';
end

% Get HK to show its log
try
    iHK = strmatchi('HK',{o.name},'exact');
    h=gr.plotMesh(log10(o(iHK).var(Iy,:,:)),'axis',ax); colormap;
    title([tts ', showing log10(HK)']);
catch %#ok
    fprintf('Can''t find HK, try TRAN\n');
    [~, HK] = o.BCF2LPF(o);
    h=gr.plotMesh(log10(HK(Iy,:,:)),'axis',ax); colormap;
    title([tts ',showing TRAN']);
end



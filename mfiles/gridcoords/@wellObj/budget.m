function [Qio,Mio] = budget(o,basename,iComp)
% B     = well.budget();          % compute budget of well over time
% B     = well.budget(basename);  % compute budget of wells over time
% [B,M] = well.budget(basename,icomp); also compute mass budget over time
%
%  TO 110822

if nargin<2 || ~exist('basename','var') || isempty(basename)
    load('name');
end

if nargin>2
    C = readMT3D(sprintf('MT3D00%d.UCN',iComp));
    %CS = readMT3D(sprintf('MT3D00%dS.UCN',iComp));  % sorbed conc
end

H=readDat([basename,'.hds']);
B=readBud([basename,'.bgt']);

time = [H.totim]; % length(B);

%% Print zone budgets for all wells

ZONE=zeros(size(B(1).term{1}));
fprintf('1) time in user units at end of each time step or stress period.\n');
fprintf('2) For each well the column of input Q (from well(iw).Q) and output Q (from budget file)\n');
fprintf('%12s','time');
for iw = numel(o):-1:1,
    ZONE(o(iw).idx) = o(iw).nr;
    for it = 1:length(time)
        Qio(2*iw,it) =sum(B(it).term{strmatchi('WELLS',B(it).label)}(ZONE==o(iw).nr));
    end
    Qio(2*iw-1,:) = o(iw).Q(:);
    fprintf(' %10s_i %10s_o',o(iw).name,o(iw).name);
end
fprintf('\n');
fmt = ['%12g' repmat(' %12g %12g',[1,numel(o)]) '\n'];
fprintf(fmt,[time; Qio]);
 Qio = Qio';

fprintf('\n\n');

%% Mass balance of requested (STILL UNDER CONSTRUCYION, WILL BE GOOD !!)
if exist('iComp','var')
    fprintf('1) time in user units at end of each time step or stress period.\n');
    fprintf('2) For each well the mass rate Q*C of input (from well(iw).C(iComp,:)) and output Q*C (from concentration file)\n');
    fprintf('%12s','time');
    for iw = numel(o):-1:1,
        for it = 1:length(time)
            Mio(2*iw,it) =sum(B(it).term{strmatchi('WELLS',B(it).label)}(ZONE==o(iw).nr).*...
                              C(it).values(ZONE==o(iw).nr));
        end
        Mio(2*iw-1,:) = o(iw).Q.*o(iw).C(iComp,:);
        fprintf(' %10s_i %10s_o',o(iw).name,o(iw).name);
    end
    fprintf('\n');
    fmt = ['%12g' repmat(' %12g %12g',[1,numel(o)]) '\n'];
    fprintf(fmt,[time; Mio]);    
end

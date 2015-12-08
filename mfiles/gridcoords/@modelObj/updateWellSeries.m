function o = updateWellSeries(o)
%% o = updateWellSeries(o) -- update well series in Model   
% Remove all wellSeries for which there are not children (wells).
%
% TO 120613

%% Get the wellSeries and the wells from te cell array o:
wellSeries  = o(strmatchi('wellSeriesObj',{o.type})).var;
well        = o(strmatchi('wellObj'      ,{o.type})).var;

%% Get the still active wellSeries (who have children among the wells)
parents = unique([well.parent]);   % Remaining well-series id's:
active = false(size(wellSeries));  % boolean list
for i=1:length(wellSeries),
    active(i)=ismember(wellSeries(i).id,parents);
end
wellSeries=wellSeries(active);  % <--cleaned well series

%% Continue with the remaining wellSeries only to update their list of children
parent=[well.parent];  % parent of each individual well
for i=1:length(wellSeries)
%    fprintf('\n');
%    fprintf('WellSeries(%2d).children =',wellSeries(i).id); fprintf('%4d',wellSeries(i).children); fprintf('\n');
    wellSeries(i).children = [well(wellSeries(i).id==parent).id];
%    fprintf('WellSeries(%2d).children =',wellSeries(i).id); fprintf('%4d',wellSeries(i).children); fprintf('\n');
end

%% Upate the o cell array for output
o(strmatchi('wellSeriesObj',{o.type})).var = wellSeries;

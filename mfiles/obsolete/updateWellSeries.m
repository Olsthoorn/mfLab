function Model = updateWellSeries(Model)
%% Model = updateWellSeries(Model)
%
% The model is a cell array with the arrays making up the model
% Each line holds one object or array
% {VariableName  Object  objectType}
% for instance
%
% {'RIV'        RIV        'stress';
%  'PEFF'       PEFF       '3Dlay';
%  'RECH'       RECH       '3Dtime':
%  'HDS'        HDS        'struct';
%  'STRTHD'     STRTHD     '3Dlay'
%  'LAYWET'     LAYWAT     'zlist';
%  'well'       well       'wellObj':
%  'wellSeries' wellSeries 'wellSeriesObj';
%  'gr'         gr         'gridObj'};
%   
% Here we remove all wellSeries for which there are not children (wells).
%
% TO 120613

%% Get the wellSeries and the wells from te cell array Model:
wellSeries  = Model{strmatchi('wellSeriesObj',Model(:,3)), 2};
well        = Model{strmatchi('wellObj'      ,Model(:,3)) ,2};

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

%% Upate the Model cell array for output
Model{strmatchi('wellSeriesObj',Model(:,3)),2}=wellSeries;

function o = subModel(o,xlim,ylim,Ilay)
% Model = Model.subModel(xlim,ylim,Ilay) -- cut submodel from a larger model
% xlim is the part of the x-axis to select (xmin xmax) of new grid.
% ylim is the part of the y-axis to select (ymin ymax) of new grid
% Ilay   are the layers of the new grid (subection of those of the old grid.
%
% TO 120607 120814

gr = o.grid();

% having the grid compute the indices of the selected model
Ix = between(gr.xm,xlim);
Iy = between(gr.ym,ylim);

if nargin<4,
    Ilay=1:gr.Nlay;
else
    Ilay = Ilay(ismember(Ilay,1:gr.Nz));
end

% Get the other arrays (backward because of memory allocation
for i=1:length(o)
    % Cutout and store the submodel for this variable
    o(i) = gr.cutout(o(i),Ix,Iy,Ilay);
end

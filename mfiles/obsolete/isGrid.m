function B = isGrid(gr)
% B = isGrid(gr) -- tests for grid
% useful with cellfunn to find objects of a class fast
% e.g.   cellfun(@isGrid,varargin)
% TO 130404

B = strcmp(class(gr),'gridObj');
function IdxMF = IdxMatlab2Modflow(IdxML,gr)
%IDXMATLAB2MODFLOW  converts global Matlab index to global MODFLOW index
%
% Example:
%    IdxMF = IdxMatlab2Modflow(IdxML,gr);
%
% SEE ALSO: cellIndices cellIndex xyzindex linegrid inpolygon

%
% TO 120823 130227

LRC = cellIndices(IdxML,gr.size,'LRC');

IdxMF = gr.Nxy*(LRC(:,1)-1)+gr.Nx*(LRC(:,2)-1)+LRC(:,3);

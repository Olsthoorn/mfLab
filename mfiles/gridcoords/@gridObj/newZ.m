function Znew = newZ(o,planes,subdivisions)
%GRIDOBJ/NEWZ --  Generate zNew, a new 3D grid from gr.Z and planes and subdivisions in gr as given
%
% USAGE:
%    Znew = newZ(o,plane,subdivisions) ---
%
% inputs:
%  planes are indices of the planes in o.Z that will be planes in the Znew
%     (except plane 1, which is always implied).
%  subdivisions are numbers indicating into how many equally thick sublayers
%     the space between the indicated planes will be subdivided.
% EXAMPLE
%   Znew = gr.newZ([2 4 7],[4,5,4]);
%      meaning that planes gr.Z(:,:,[1 2 4 7]) will be in newZ while their
%      intermdiate space is subdibided into 4, 5 and 4 equally thick layers
%      repectively.
% TO 120830 

if planes(end) ~= o.Nlay+1
    error('%s: the last plane Nr (%d) must equal the last plane in the grid = %d',mfilename,planes(end),gr.Nz+1);
end
if planes(1) == 1
    error('%s: the first plane must be >1, as 1 is always implied',mfilename);
end
if numel(planes)~=numel(subdivisions)
    error('%s: Number of planes (%d) must equal number of subdivisions (%d)',mfilename,numel(planes),numel(subdivisions));
end
if any(planes>o.Nz+1) || any(planes<1)
    error('%s: planes must be increasing and all smaller than the number of layers in the olde grid %d',mfilename,o.Nz+1);
end
if numel(planes)>1 && ~all(diff(planes))>0
    error('%s: repeated planes or a non ascending plane-snumber sequence is not allowd',mfilename);
end
if any(subdivisions<1),
    error('%s: all subdivisions must be >=1',mfilename);
end

planes               = [1 planes(:)'];
Nnew                 = sum(subdivisions);
Znew                 = NaN(o.Ny,o.Nx,Nnew);
Znew(:,:,1)          = o.Z(:,:,planes(1)); % planes(1) is always 1

k=1;
for ip=2:length(planes)
    dz = max(o.MINDZ,(o.Z(:,:,planes(ip-1))-o.Z(:,:,planes(ip)))/subdivisions(ip-1));
    for j=1:subdivisions(ip-1)
        Znew(:,:,k+1)=Znew(:,:,k) - dz;
        k=k+1;
    end
end

fprintf \n


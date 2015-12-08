function showStress(o,stressName,stress,ttl,iPer,column,colHdr)
% gridObj/showStress(stressName,stress[,ttl[,iPer[,col[,colHdr]]]) -- contours stresses
%
% stressName: 'RIV', 'GHB', 'DRN', 'STR', 'DRT', 'PNTSRC' 'WEL'
% stress: a list or array if cells each with a list of the form
%    of the stresses RIV GHB and, i.e. with [iPER,LAY,ROW,COL,variables] on
%    each line.
% ttl:     start of title on graph (often just use basename in the call)
% iPer:   stress period number
% column: column of list to be shown
% colHdr: title of column of list to show in title
%   stress
%   You can specify the column, if omitted, col 5 is used (usually the head)
%   if iPer (stress period number) is omitted, 1 is used.

if nargin<6, column = 5;    end
if nargin<5, iPer=1;        end
if nargin<4, ttl='???';     end

if nargin<7
    switch stressName
        case 'RIV'
            hdr = {'stress period','layer','row','col','river stage','river cond','river bottom elev.'};
        case 'GHB'
            hdr = {'stress period','layer','row','col','prescr. head','GHB cond'};
        case 'DRN'
            hdr = {'stress period','layer','row','col','drain elevation','drain conductance'};
        case 'WEL'
            hdr = {'stress period','layer','row','col','presribed flow'};
        case 'PNTSRC'
            hdr = {'stress period','layer','row','col','C1','ITYPE','C1','C2','C3','C4'};
        otherwise
            hdr = {'stress period','layer','row','col','column5','column6','column7','column8'};
    end
    colHdr = hdr{max(1,min(length(hdr),column))};
end

if iscell(stress), stress = cell2list(stress); end

if size(stress,2)<5, error('%s: This can''t be a legal stress like RIV, GHB etc',mfilename); end

stress = stress(stress(:,1)==iPer,:);

var      = o.const(NaN);
Idx      = cellIndex(stress(:,[4,3,2]),o.size);
var(Idx) = stress(:,column);

range = ContourRange(var,50);
caxis([min(range),max(range)]);

LAY    = unique(stress(:,2));

for iLay = LAY(:)'
    figure; hold on; xlabel('x [m]'); ylabel('y[m]'); zlabel('z[m]');
    title(sprintf('%s: %s, %s for stress period %d and layer %d',ttl,stressName,colHdr,iPer,iLay));

    h(iLay) = surf(o.XGr,o.YGr,o.ZGrBlay(:,:,iLay),var(:,:,iLay),'edgecolor','none');
    
    colorbar;
end

function mf_lay2mdl(basename,xGr,yGr,varargin)
%MF_LAY2MDL -- gets LAYER info from SHEET LAY in accompanying workbook
% and puts the corresponding 3D layer arrays in the workspace.
% The function simplifies geneation of a generic model.
%
% It only generates the arrays
%   HK, VK, SY, SS, PEFF, VKCB, PORCB and the grid object
%   
%   generate SY and SS only when transient
%   as deduced from the column transient in sheet PER
%
%   generate VKCB and PORCB only when there is a confining bed (any (LAYCBD) ~=0).
%
% USAGE:
%   mf_lay2mdl(basename,xGr,yGr,'Nz',12,varargin)
%   mf_lay2mdl(basename,xGr,yGr,'Nz',20,'aL','DMCOEF_1','DMCOEF_1','RHOB');
%
%  varargin may contain any headers in the LAY worksheet between braces as in
%  usage example 2, where each name generates a 3D model array of that
%  name.
%
% TO 130612

% LAYER	Top/d	LAYCBD	HK	VK	Sy	Ss	Peff strthd'

[LAYnams,LAYvals] = getLayers(basename);
[PERnams,PERvals] = getExcelData(basename,'PER','H');

transient = any(PERvals(:,strmatchi('transient',PERnams)));

%% Always needed

layInfo =  LAYvals(:,strmatchi('LAYCBD',LAYnams));
LAYCBD  =  layInfo;
LAYCBD(find(layInfo)-1) = 1;
LAYCBD(     layInfo==1) =[];

clear layInfo

zGr = LAYvals(:,strmatchi('D',LAYnams,'exact'));
zGr = [0; -cumsum(zGr)];  % 0 = defaultTop of model, change it in mf_adapt

[zTop,varargin] = getProp(varargin,'zTop',[]);
if isempty(zTop)
    i = strmatchi('top',LAYnams,'exact');
    if i(1)
      zTop = LAYvals(1,i(1));
    else
        zTop = 0;
    end
end

gr = gridObj(xGr,yGr,zGr+zTop,'LAYCBD',LAYCBD);

%% Get basic 3D arrays that are generally needed for MODFLOW
HK      = LAYvals(gr.ITlay,strmatchi('HK'  ,LAYnams,'exact')); %#ok
VK      = LAYvals(gr.ITlay,strmatchi('VK'  ,LAYnams,'exact')); %#ok
STRTHD  = LAYvals(gr.ITlay,strmatchi('STRTHD',LAYnams,'exact')); %#ok

if transient
    SS      = LAYvals(gr.ITlay,strmatchi('SS'  ,LAYnams,'exact')); %#ok
    SY      = LAYvals(gr.ITlay,strmatchi('SY'  ,LAYnams,'exact')); %#ok
end

PEFF    = LAYvals(gr.ITlay,strmatchi('PEFF',LAYnams,'exact')); %#ok

if any(gr.LAYCBD)
    VKCB    = LAYvals(gr.ITcbd,strmatchi('VK'  ,LAYnams,'exact')); %#ok
    PORCB   = LAYvals(gr.ITcbd,strmatchi('PEFF',LAYnams,'exact')); %#ok
end


%% Get other requested arrays
for i=1:numel(varargin)
    j = strmatchi(varargin{i},LAYnams);
    if ~j(1)
        error('%s: Requested header <<%s>> does not exist in sheet LAY of workbook <<%s>>',...
            mfilename,basename);
    end
    if length(j)>1
        error('%s: Requested header <<%s>> is not unique in sheet LAY of workbook <<%s>>',...
            mfilename,basename);
    end
    eval([varargin{i} '= LAYvals(:,strmatchi(varargin{i},LAYnams));']); 
end

clear basename i j transient varargin LAYnams LAYvals PERnams PERvals xGr yGr zGr layInfo LAYCBD


vars = setdiff(who(),{'gr','zTop'});

for i = 1:numel(vars)
    if ~exist(vars{i},'var') || isempty(eval(vars{i}))
        error(['%s: missing variable <<%s>>,mfilename.\n',...
             'REMEDY: add column <<%s>> to sheet LAY'],mfilename,vars{i},vars{i});
    end
    eval([vars{i} '= gr.const(' vars{i} ');']);
end

IBOUND  = gr.const(1); %#ok

clear i vars

% Show what is save to the workspace
fprintf('%s: The following arrays and variables are saved to the Matlab workspace:\n',mfilename);
w = whos();
fprintf('%20s\t%s\t%9s\t%s\n','name','size','class','bytes');
for i=1:numel(w)
    fprintf('%20s',w(i).name);
    fprintf( '\t%d %d %d',w(i).size);
    fprintf('\t%10s',w(i).class);
    fprintf('\t%d\n',w(i).bytes);
end
%display(whos());


%% save local variables to workspace
save2base(true);
function [bcn,varname]= getBCN(o,varname,BCN)
%% [bcn,varname]= gr.getBCN(o,type,var) -- 3D array with boundary conditions
% i.e. from 'RIV', 'GHB' or 'DRN'
% BCN is the corresponding variable RIV, GHB or DRN
%
% type can be [L][C|H]]RIV|GHB|DRN
% where means 10log of varible, C means conductance, H means head
%
% for use in showLayers, where one could fill into the edit box
%   LCriv Hriv Cdrn HGHB hGHB Hghb L[H]ghb (head is default
%   so the string must have 3, 4 of 5 characters to work.
% The function assumes that the RIV, DRN and or GHB can be loaded
% into the workspace throug load(basename). This loading is done
% internally in showLayers.
%
% TO 120604

if iscell(BCN), BCN = cell2list(BCN); end

if nargin<2|| ~ischar(varname) || length(varname)<3
    warning('gridObj:getBCN:insufficientOrWrongInputArguments',...
        '%s: insufficient or wrong arguments, see help of this mfile',mfilename);
    bcn=[];
    return;
end

varname=upper(varname);

%% See if logarithm is desired
if length(varname)>4
    logscale = upper(varname(1))=='L';
    varname=varname(2:end);
else
    logscale=0;
end

%% See if conductance is desired
if length(varname)>3
    if varname(1)=='C',
        varname= ['conductance of ' varname(2:end)];
        column = 6;  % conductance
    else % assume head whatever user inputs
        varname= ['head of ' varname(2:end)];
        column =5;   % head
    end
    
else
    column=5;  % head (default)
end

if logscale, varname =['logarithm of ' varname]; end

%%
% Global array index of the specified boundary condition (stress)
% Only look at values for the first stress period
Ibcn = cellIndex(BCN(BCN(:,1)==1,[4,3,2]),o.size);

%% 3D array to show this stress
bcn= o.const(0);

if logscale
    % Use only first stress period
    bcn( Ibcn) = log10(BCN(BCN(:,1)==1,column));
else
    % Use only first stress period
    bcn( Ibcn) =       BCN(BCN(:,1)==1,column) ;
end

% Unique layers in which this BCN
layersWithBCN = unique(BCN(:,2));
for iLay=1:size(bcn,3)
    if ~ismember(iLay,layersWithBCN)
        bcn(:,:,iLay)=NaN;
    end
end
%% Start varname expression with capital letter
varname(1)=upper(varname(1));


function [LS,c,m,L] = isLineSpec(lineSpec)
%ISLINESPEC check if lineSpec is a legal Matlab lineSpec
%
% Example:
%    [LS,c,l,m] = isLineSpec(lineSpec) -- 
%
%    legal lineSpecs are a combination of
%    a color from       'bwgkmcyw'
%    a line type from   '--',':','-.','-','none'
%    a marker from      '.os^vp+x*'
%
% see also: mf_color mf_lineType isColor
%
% TO 011021 151203

LSPECS  = {'--',':','-.','-'};
COLORS  = 'brgkmcyw';
MARKERS = '.os^vp+x*';

LS = false; c  =''; m  =''; L  ='';

if ~ischar(lineSpec)
    if isa(lineSpec,'double') && numel(lineSpec)==3
        c = lineSpec(:)';
        LS = true;
    end
    return;
end

if strcmpi(lineSpec,'none')
    c = lineSpec;
    LS = true;
    return;
end

%% Does lineSpec have colors ?
Lclr      = ismember(lineSpec,COLORS);
c        = lineSpec( Lclr(1));

%% Does lineSpec have markers ?
Lmrk      = ismember(lineSpec,MARKERS);
m        = lineSpec( Lmrk(1));

%% Does lineSpec have lineStyles ?
for i=1:numel(LSPECS)
    j=strfind(lineSpec,LSPECS{i});
    if ~isempty(j)
        L = LSPECS{i};
        break;
    end
end

LS = ~isempty(c) || ~isempty(m) || ~isempty(L);
    

function [LC,clr] = isColor(varargin)
%ISCOLOR check if colorSpec is a legal Matlab color
%
% Example:
%    [LC, clr] = isColor(color) -- 
%
%     LC is true or false
%     clr is the color if LC is true else is []
%
%    legal colors are one of 'brwgkmcyw'
%    or a RGB triple   [R G B; R G V; ...]
%
% see also: mf_color mf_colorType isLineSpec
%
% TO 130825

[clr,varargin] = getNext(varargin,'char',[]);
[RGB,       ~] = getNext(varargin,'double',[]);

if ~isempty(clr)
    LC = all(ismember(clr,'brgkmcyw'));
    if LC, clr = clr(:); else clr=[]; end
    return;
elseif ~isempty(RGB)
        LC = size(RGB,2)==3 && all(RGB(:)>=0 & RGB(:)<=1);
        if LC, clr = RGB; else clr=[]; end
else
    LC = false;
    clr = [];
end
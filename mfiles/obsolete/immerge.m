function out = immerge(bg, fg, coef)
% Creates an image 'out' of same type as bg.
% 'out' has 'bg' as background, and 'fg' (transparent, 
% weighted by 'coef') above 'bg'.
% 'out', 'bg', and 'fg' are RGB images.
%
% merged = immerged(bg, fg, coef)
%	- bg matrix of type double or uint8
%	- fg matrix of type double or uint8
%	- coef is a scalar between 0 and 1, or a matrix of 
%	  such scalars, same size as 'fg' and 'bg' (AlphaData).
%
% Gauthier Fleutot 28-07-2004
% fleutotg@esiee.fr

if ~all(size(coef) == [1 1])
	coef = repmat(coef,[1 1 3]);	% extend the coef matric in the 3rd dim.
end
out = coef.*bg + (1-coef).*fg;

%out = uint8(out);

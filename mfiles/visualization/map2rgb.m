function RGB = map2rgb(X,map)
%MAP2RGB converts flat image with colormap map to RGB
%
% Example:
%    RGB = map2rgb(X,map); 
%
% TO 121017

R = map(:,1);
G = map(:,2);
B = map(:,3);

RGB = cat(3,R(X+1),G(X+1),B(X+1));

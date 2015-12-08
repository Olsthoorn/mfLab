function clayers(xm,ym,VAR,varname)
%CLAYERS contours all layers where VAR is the 3D-value matrix
%
% Example:
%    CLAYERS(xm,ym,VAR,varname)
%
% TO 090724

figure;
for i=1:size(VAR,3)
    subplot(2,2,i)
    contour(xm,ym,VAR(:,:,i))
    fprintf('Layer %d: %smin=%12g, %smax=%12g\n',i,...
        varname,min(min(VAR(:,:,i))),...
        varname,max(max(VAR(:,:,i))));
    title(sprintf('%s, layer %d',varname,i));
end

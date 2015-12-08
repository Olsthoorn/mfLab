function map=mf_colormap(clrs,L)
%MF_COLORMAP make your own colormap of size [L,3] using colors 'b r g k m c y w'
%
% Example:
%    map=mf_colormap(clrs,L)
%    map=mf_colormap('mbry',64);
%
%  make your own colormap of size [L,3] using colors 'b r g k m c y w'
%  specified in a string like 'bcyg'  'bryw' etc
%
% See also: mf_color
%
% TO 120409

for ic=1:length(clrs)
    switch clrs(ic)
        case 'w',  c=[1 1 1];
        case 'b',  c=[0 0 1];
        case 'c',  c=[0 1 1];
        case 'm',  c=[1 0 1];
        case 'y',  c=[1 1 0];
        case 'r',  c=[1 0 0];
        case 'g',  c=[0 1 0];
        case 'k',  c=[0 0 0];
        otherwise % white
            c=[1 1 1];
    end
    C(ic,:)=c;
end

ix=1:L;
n=length(clrs);
dL=(L-1)/(n-1);
x=1+(0:(n-1))*dL;

map=interp1(x,C,ix);
function o = round(o,sz)
% well = wel.round(sz) -- rounds x and y coordinates of well to sz
%
% TO 120826

for iw = 1:length(o)
    o(iw).x = round(o(iw).x/sz)*sz;
    o(iw).y = round(o(iw).y/sz)*sz;
end
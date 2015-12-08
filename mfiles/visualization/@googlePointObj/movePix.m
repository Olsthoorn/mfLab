function o = movePix(o,npx,npy)
    % o = o.movePix(npx,npy);
    % move point by npx,npy pixels using zoom level op point
    % setUR and LL of object in LL from center, zoom and size
    
    [Lat Lon] = pix2LL(o.ix,o.iy,o.px+npx,o.py+npy,o.zoom);
    
    o = googlePointObj(Lat,Lon,o.zoom);
    
end

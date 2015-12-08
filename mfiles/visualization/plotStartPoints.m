function h = plotStartPoints(fname,lSpec)
%PLOTSTARTPOINTS plots the endpoints in MODPATH6 endpoint file named fname
%
% Example:
%    hdl = plotStartPoints(fname,[lSpec])
%    where lSpec is a lineSpec as valid for plot(x,y,lSpec)
%
% See also: plotEndPoints plotTsrPoints plotPath
%
% TO 130219
    
    endp = readEndPoints(fname);
    
    ix = strmatchi('xGI',endp.colHdr);
    iy = strmatchi('yGI',endp.colHdr);
    iz = strmatchi('zGI',endp.colHdr);
    h(endp.groupCount) = 0;
    for iGrp = 1:endp.groupCount
        I = find(endp.P(:,2)==iGrp);
        if nargin<2
            lSpec = [mf_color(iGrp,'brgkmc'),mf_marker(iGrp,'o*^sp+x')];
        end
        h(iGrp) = plot3(endp.P(I,ix),endp.P(I,iy),endp.P(I,iz),lSpec);
    end
end

function h = plotEndPoints(fname,lSpec)
%PLOTENDPOINTS plots the endpoints in MODPATH6 endpoint file named fname
%
% Examle:
%    hdl = plotEndpoints(fname,[lSpec])
%    where lSpec is a lineSpec as valid for plot(x,y,lSpec)
% 
% See also: plotPath plotTsrPoints plotStartPoints
%
% TO 130219
    
    endp = readEndPoints(fname);
    
    ix = strmatchi('xGF',endp.colHdr);
    iy = strmatchi('yGF',endp.colHdr);
    iz = strmatchi('zGF',endp.colHdr);
    h(endp.groupCount) = 0;
    for iGrp = 1:endp.groupCount
        I = find(endp.P(:,2)==iGrp);
        if nargin<2
            lSpec = [mf_color(iGrp,'brgkmc'),mf_marker(iGrp,'o*^sp+x')];
        end
        h(iGrp) = plot3(endp.P(I,ix),endp.P(I,iy),endp.P(I,iz),lSpec);
    end
end

function hdl = plotPath(pathFileName,~)
%PLOTPATH  plots the MODPATH6 path lines of the particle groups
%
% Example:
%     hdl = mpath_particleGroupObj.plotPthPoints(pathFileName[,lSpec])
%     where lSpec is a lineSpec as valid for plot(x,y,lSpec)
%
% See also: plotEndPoints plotStartPoints plotTsrPoints
%
% TO 130219
    
    pth = readPath(pathFileName);
  
    % sorting according to group,particleID,cumulativeTimeStep
    pth.P = sortrows(pth.P,[2,1,4]);
    
    ix = strmatchi('xG',pth.colHdr);
    iy = strmatchi('yG',pth.colHdr);
    iz = strmatchi('zG',pth.colHdr);

    % now ordered
    [~,m1] = unique(pth.P(:,1),'first');
    [~,m2] = unique(pth.P(:,1),'last');

    hdl(numel(m1),1) = NaN;
    for iP = 1:numel(m1)
        Irange = m1(iP):m2(iP);
        lSpec   = mf_color(pth.P(m1(1),2),'rbgkmcy');
        hdl(iP) = plot3(pth.P(Irange,ix),...
                        pth.P(Irange,iy),...
                        pth.P(Irange,iz),...
                                lSpec);
    end

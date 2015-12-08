%% determineGradGroundSurface.m ---
% determine the surface gradient in a given area
% desired to set the gradient in the fixed-gradient boundary option.
%
% TO 141107
clear xMid yMid dX dY grX grY

% get the dX, dY and dZ to compute the gradients.
xMid = 0.5* (gr.Xm(:,1:end-1) + gr.Xm(:,2:end)); dX = diff(gr.Xm,1,2); dZx = diff(gr.Z(:,:,1),1,2);
yMid = 0.5* (gr.Ym(1:end-1,:) + gr.Ym(2:end,:)); dY = diff(gr.Ym,1,1); dZy = diff(gr.Z(:,:,1),1,1);

% ground-surface gradient in both directions
grX = (dZx./dX); grX_ = 0.5*(grX(1:end-1,:)+grX(2:end,:));
grY = (dZy./dY); grY_ = 0.5*(grY(:,1:end-1)+grY(:,2:end));

% total gradient
gr_  = sqrt(grX_.^2+grY_.^2);

% locations of gradient
xM_ = 0.5*(xMid(1:end-1,:)+xMid(2:end,:));
yM_ = 0.5*(yMid(:,1:end-1)+yMid(:,2:end));

% show gradient map, selected areas
figure; hold on; contour(xM_,yM_,gr_,50);
plot(gradAreaNorth(:,1),gradAreaNorth(:,2),'r','lineWidth',2);
plot(gradAreaSouth(:,1),gradAreaSouth(:,2),'r','lineWidth',2);

% get points within the respective north and south areas
INnorth = inpolygon(xM_,yM_,gradAreaNorth(:,1),gradAreaNorth(:,2));
INsouth = inpolygon(xM_,yM_,gradAreaSouth(:,1),gradAreaSouth(:,2));

% compute the gradient in these areas
grNorth = median(gr_(INnorth));
grSouth = median(gr_(INsouth));

% and print them
fprintf('gradient ground surface north = %.4g\n',grNorth);
fprintf('gradient ground surface south = %.4g\n',grSouth);

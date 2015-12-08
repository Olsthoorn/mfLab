function [xGr,yGr] = wgs2utmContinuous(varargin)
   %WGS2UTM_CONTINUOUS - get xGr and yGr of a grid defined by Egr,NGr but make
   % sure that xGr is continuous across vertical UTM 6 deg boundaries that separate
   % UTM zones.
   %
   % USAGE: [x,y] = wgs2utmContinuous(N,E[,Nrange,Erange,][clr])
   %
   % N and E input can define a grid or a set of points. N and E are
   % considered a set of points if ~all(size(N)==size(E)).
   % To make sure that it works over desired coordinates, use Nrange and
   % Erange to give the total span of the desired coordinates. This is necessary
   % because individual calls have no idea about the span over all
   % coordinates. If you omit Nrange or Erange the N and E are used to
   % figure out the span, but this may not suffice your needs.
   %
   % Explanation: The UTM zones are 6 degrees wide. If your grid spans a
   % wider zone than 6 degrees or crosses a 6 degree boundary, you would
   % end up with discontinuous x-coordinates. To prevent that, the
   % coordinates will be made continuous by adding the end x-coordinate of
   % the previous zone to the next etc. until in the right-most zone.
   % The x-coordinates of the left-most zone are the true UTM coordinates
   % for that zone.
   % The NGr range is necessary to compromize the connection between
   % adjacent zone with the least distortion possible.
   %
   % TO 131212
   
   [Ngr   ,varargin] = getNext(varargin,'double',[]);
   [Egr   ,varargin] = getNext(varargin,'double',[]);
   [Nrange,varargin] = getNext(varargin,'double',[min(Ngr) max(Ngr)]);
   [Erange,varargin] = getNext(varargin,'double',[min(Egr) max(Egr)]);
   [clr   ,  ~     ] = getNext(varargin,'char',[]);

   if ~isempty(clr), hold on; end
   
   if all(size(Egr)==size(Ngr)) % assume coordinate pairs
         [  ~,yGr] = wgs2utm(Ngr,Egr);
   else  % assume coordinates of a grid
         [  ~,yGr] = wgs2utm(Ngr,mean(Egr(:))*ones(size(Ngr)));
   end

   Ngr=Ngr(:)'; Egr=Egr(:)';
  
   UTM_Eboundaries = (0:6:360) - 180;
   UTM_Eboundaries = UTM_Eboundaries(UTM_Eboundaries<=max(Erange) & UTM_Eboundaries>= min(Erange));

   if isempty(UTM_Eboundaries)
       [xGr,  ~] = wgs2utm(mean(Nrange)*ones(size(Egr)),Egr);
       if ~isempty(clr)
           plot(xGr,yGr,clr);
       end
       return;
   end
   % use mean northing of grid
   N_Avg = mean(Nrange(:));

   % get the UTM coordinates for these boundaries at mean northing
   % but subtract small value to get highest x-coordinate at this
   % northing within UTM band
   [xGrBoundaries_min, ~] = wgs2utm(N_Avg*ones(size(UTM_Eboundaries)),UTM_Eboundaries+1e-8);
   [xGrBoundaries_max, ~] = wgs2utm(N_Avg*ones(size(UTM_Eboundaries)),UTM_Eboundaries-1e-8);
    xGrBoundaries_max     = cumsum(xGrBoundaries_max);
   [xGr          , ~]     = wgs2utm(N_Avg*ones(size(Egr))            ,Egr);

   % recursively increase xGr of successive blocks with highest UTM values of
   % blocks to their left, to obtain continuously increasing x-coordinates

   for i=1:numel(UTM_Eboundaries)
       xGr(Egr>UTM_Eboundaries(i)) = xGr(Egr>UTM_Eboundaries(i))+(xGrBoundaries_max(i)-xGrBoundaries_min(i));
   end
   
   if ~isempty(clr)
       plot(xGr,yGr,clr)
   end   
end

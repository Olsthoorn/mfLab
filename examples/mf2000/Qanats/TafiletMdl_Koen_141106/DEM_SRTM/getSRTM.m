%% getSRTM -- extracts a set of subdems of diffrent scale from a 3arc second (90 m) SRTM dem
%
% for the current Tafilalet model we use dem 
%
%     srtm_36_06.asc
%
% The reslution factor can be 1, 2, 3, etc .... meaning that the resolution
% of the output DEM is based on 1x1, 2x2, 3x3 etc points of the original
% DEM, its scale wil thus be 90, 180, 270, 360 m or about. Actually 3, 6,
% 9, 12 etc arc seconds.
% The resolution will lbe in the name of the generated DEM file expressed
% in the approximate cell size as written above.
%
% The resulting DEM will also be smoothed using medfilt
%
% The resulting files will be mat file for efficience, that contain
% The DEM and its coordinatea in Lat Lon and X, Y (utm for the arc range)
%
% TO 131010

%% Background info.
%  See readme.txt for background info from the producer.
%  Interpretation. One STRM file is 5 deg by 5 deg subdivided into 6001 by
%  6001 points. Notice from the header
%
% ncols         6001
% nrows         6001
% xllcorner     -5.0004160855315
% yllcorner     29.999583309109
% cellsize      0.00083333333333333
% NODATA_value  -9999
%
% That the left coordinate is half the cell size left of -5deg
% That the lower coordinates is half the cell size below 30deg.
% Hence, the points are exactly on 5 and 0 deg and 30 and 35 deg.
% The itnermediate distance is the cellsie. 6000 cells fit between the
% boundaries of the STRM block. The total number of points is indeed 6001
% by 60001. The boundary points should exactly coincide with the boundary
% points of adjacent blocks.
%
% The size of one cell is 5*60*60/(6000)=3 seconds
% The vertical distance of once cell is then approximately
% 3/(60*60*360)*40000000 m =  92.6 m
% The same horizontally at about 31.5 deg is
% 92.6 * cos(31.5/180) = 91.2 m.
%
% The actual values are integer and therefore in meters. The accuracy may
% not be very high depending on the flatness of the terrain. It may vary
% between 1 m for lauchnpads and salt lakes to 16 m in accidented terrain.
%
% The actual elevaton of a country may also differ from the WGS84
% elevation.
%
% See http://digaohm.semar.gob.mx/imagenes/hidrografia/S60_Ed3Eng.pdf:
% The most accurate approach for obtaining WGS 84 coordinates is to acquire satellite tracking data
% at the site of interest and position it directly in WGS 84 using GPS positioning techniques. Direct
% occupation of the site is not always possible or warranted. In these cases, a datum transformation
% can be used to convert coordinates from the local system to WGS 84.
%
% Morocco uses the Merchick local geodetic system with Clark 1880
% associated ellipsoid (code CD).
% Transformation parameters are DX = 31+/-5m, DY=147+/-3m, DZ=47+/-3m.
% As WGS is defined as WGS = local + correction
% This implies that we have to add 47 m to local z-coordinates to match
% them with WGS84. We can neglect the horizontal corrections in our case.
%
% For me it is still unclear if and how these values should be applied.
%
%

%% Domain boundaries for which a DEM is to be cut out at full resolution
E = [-4.680 -4.300];
N = [31.365 31.600];

%% Get the original DEM

FName = 'srtm_36_06.asc';

fid = fopen(FName);
for i=1:6
    s = fgetl(fid); iSpace = find(s==' '); iSpace=iSpace(1);
    s(iSpace(1)) = '='; s=[s ';']; %#ok
    eval(s);
end

Data = fscanf(fid,'%d',[ncols,Inf])';

fclose(fid);

%% check the number of rows obtained
if size(Data,1)~=nrows
    error('nrows = %d expected, nrows = %d read from file %s',...
        nrows,size(Data,1),FName);
end


%% Extract the model area

% The coordinates of any pixcel center can be found from
eastGr = xllcorner + (0:ncols)*cellsize;
northGr= yllcorner + (0:nrows)*cellsize;

% Important: Because the raster data are such that the first line is the top
% line in the raster, i.e. the highess y coordinate, we must sort the
% northC coordinates so that they are descending
northGr = sort(northGr,'descend');

%% Cell centers
eastC = 0.5*(eastGr( 1:end-1)+eastGr( 2:end));
northC= 0.5*(northGr(1:end-1)+northGr(2:end));

%% get the subset
Ix = eastC>= E(1)-0.5*cellsize & eastC <=E(end)+0.5*cellsize;
Iy = northC>=N(1)-0.5*cellsize & northC<=N(end)+0.5*cellsize; 

eastC = eastC( Ix); xllcornerNew = min(eastC )-cellsize/2;
northC= northC(Iy); yllcornerNew = min(northC)-cellsize/2;

%% Here is the new subDEM
DataNew = Data(Iy,Ix);

% For selecton of subDEMs turn DEM and northC upside down to keep
% the LL corner upon selection of subdem elevatons.
% Turn back when ready.
DataNew    = DataNew(end:-1:1,:);  % turn upside down to fix the LL corner coordinates
northC     = northC(end:-1:1);

% Issue a DEM of the desired scale resolutionFactor, which can be 1,2,3,4 what ever
% hte number of original DEM cells to combine in a new DEM where n means
% that n*n cells will be combined. 1 is the original DEM cutout on the
% desired wgs  coordinates


for resolutionFactor=1:4

    cellWidth = cellsize*resolutionFactor/2;
    
    % Select DEM but keep the LL point fixed
    DEM = DataNew(1:resolutionFactor:end,1:resolutionFactor:end);
    E = eastC( 1:resolutionFactor:end);
    N = northC(1:resolutionFactor:end);

    DEM = DEM(end:-1:1,:);
    N   = N(end  :-1:1);
   
    % the cell size varies over the grid, and so the y lines are not completly
    % parallel in the plane. As a compromize we take X and Y for the center of
    % our subgrid.
    EGr = [E(1)-cellWidth/2 0.5*(E(1:end-1)+E(2:end)) E(end)+cellWidth/2];
    NGr = [N(1)+cellWidth/2 0.5*(N(1:end-1)+N(2:end)) N(end)-cellWidth/2];
    
    [xGr,~] = wgs2utm(mean(NGr)*ones(size(EGr)),EGr);
    [~,yGr] = wgs2utm(NGr,mean(EGr)*ones(size(NGr)));
    
    
    %% Smoothing the DEM
    spans=[7 5 3 3]; degree=1;
    for k=1:6
        span = spans(resolutionFactor);
        % We first smooth in the y-direction and then in the x-direction.
        % For the y-direction:
        % Take the original DEM (Layers{1}), turn the alternating columns up side
        % down and make it a single long vector. This ensures that the data in this
        % vector will be continuous at the bottom and top of the DEM. Then apply
        % smooth on this very long vector to get smoothed values. When done reshape
        % the DEM back, and turn the altenating columns back into their original order.
        DEM1 = DEM;
        DEM1(:,2:2:end) = DEM1(end:-1:1,2:2:end);
        DEM1 = reshape(smooth(DEM1(:),span,'sgolay',degree),size(DEM));
        %DEM1 = reshape(smooth(DEM1(:),span,'loess'),size(Layers{1}));
        DEM1(:,2:2:end) = DEM1(end:-1:1,2:2:end); % restore DEM1

        % For the x-direction, first transpose the DEM, and then do the same thing
        % as before. When done back-transpose.
        DEM2 = DEM';
        DEM2(:,2:2:end) = DEM2(end:-1:1,2:2:end);
        DEM2 = reshape(smooth(DEM2(:),span,'sgolay',degree),size(DEM'));
        DEM2(:,2:2:end) = DEM2(end:-1:1,2:2:end); % restore DEM2
        DEM2 = DEM2';

        % Finally take the average of the y-smoothed and x-smoothed DEM as the
        % ultimate dem.
        DEM = 0.5*(DEM1+DEM2);

    end
    
   % DEM = medfilt2(DEM,[5,5],'symmetric');  % smooth DEM a little more

    fname = sprintf('strm_36_06_Jorf_%03dm.mat',resolutionFactor*90);
    save(fname,'DEM','EGr','NGr','xGr','yGr');
    
end

return

%% TEST

clear DEM E N X Y
for resolutionFactor=1:4
    
    fName = sprintf('strm_36_06_Jorf_%03dm.mat',90*resolutionFactor);
    
    load(fName);

    figure('name',fName,'pos',screenPos(0.75));
    axes('nextplot','add');
    xlabel('xUTM [m]'); ylabel('yUTM [m]');
    title(fName,'interpreter','none');
    
    hrange = contRange(DEM,50);    
    
    xm = 0.5*(xGr(1:end-1)+xGr(2:end));
    ym = 0.5*(yGr(1:end-1)+yGr(2:end));
    xm([1 end]) = xGr([1 end]);
    ym([1 end]) = yGr([1 end]); 
    
    [c,h] = contourf(xm,ym,DEM,hrange);
    hb = colorbar; set(get(hb,'title'),'string','elevation');
    
end


%% Save DEM model area 90 m resolutie as ESRI ascii file
% 
% FName = 'strm_36_06_Jorf_090m.asc';
% 
% fp = fopen(FName,'w');
% fprintf(fp,'ncols         %d\n'   ,size(DataNew,2));
% fprintf(fp,'nrows         %d\n'   ,size(DataNew,1));
% fprintf(fp,'xllcorner     %.13f\n',xllcornerNew);
% fprintf(fp,'yllcorner     %.13f\n',yllcornerNew);
% fprintf(fp,'cellsize      %.17f\n',cellsize);
% fprintf(fp,'NODATA_value  -9999\n');
% for iR = 1:size(DataNew,1)
%     fprintf(fp,'% d',DataNew(iR,:)); fprintf(fp,'\n');
% end
% fclose(fp);

%% smoothing the new DEM

%DEM = medfilt2(DEM,[5,5],'symmetric');  % smooth DEM a little more

%% Save DEM model area 180 m resolutie as ESTI ascii file
% 
% 
% FName = 'strm_36_06_Jorf_180m.asc';
% 
% cellsize2     = cellsize*2;
% eastC2        = eastC(1:2:end);
% northC2       = northC(1:2:end);
% DataNew2      = DataNew(1:2:end,1:2:end);
% xllcornerNew2 = eastC2(   1)-cellsize2/2;
% yllcornerNew2 = northC2(end)-cellsize2/2;
% 
% fp = fopen(FName,'w');
% fprintf(fp,'ncols         %d\n',size(DataNew2,2));
% fprintf(fp,'nrows         %d\n',size(DataNew2,1));
% fprintf(fp,'xllcorner     %.13f\n',xllcornerNew2);
% fprintf(fp,'yllcorner     %.13f\n',yllcornerNew2);
% fprintf(fp,'cellsize      %.17f\n',cellsize2);
% fprintf(fp,'NODATA_value  -9999\n');
% for iR = 1:size(DataNew2,1)
%     fprintf(fp,'% d',DataNew2(iR,:)); fprintf(fp,'\n');
% end
% fclose(fp);
% 
%% Save DEM model area 360 m resolution as ESRI ascii file
% 
% FName = 'strm_36_06_Jorf_360m.asc';
% 
% cellsize4     = cellsize*4;
% eastC4        = eastC(1:4:end);
% northC4       = northC(1:4:end);
% DataNew4      = DataNew(1:4:end,1:4:end);
% xllcornerNew4 = eastC2(   1)-cellsize4/2;
% yllcornerNew4 = northC2(end)-cellsize4/2;
% 
% fp = fopen(FName,'w');
% 
% fprintf(fp,'ncols         %d\n',size(DataNew4,2));
% fprintf(fp,'nrows         %d\n',size(DataNew4,1));
% fprintf(fp,'xllcorner     %.13f\n',xllcornerNew4);
% fprintf(fp,'yllcorner     %.13f\n',yllcornerNew4);
% fprintf(fp,'cellsize      %.17f\n',cellsize4);
% fprintf(fp,'NODATA_value  -9999\n');
% for iR = 1:size(DataNew4,1)
%     fprintf(fp,'% d',DataNew4(iR,:)); fprintf(fp,'\n');
% end
% fclose(fp);


% this yields the requested SRTM raster which just encompasses our model
% area

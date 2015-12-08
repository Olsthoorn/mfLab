function plotshp(shape)
%PLOTSHP plots the shape
%
% USAGE:
%    plotshp(shape)
%
%    plot shape created with readshp
%    'http://www.esri.com/library/whitepapers/pdfs/shapefile.pdf'
%
%    shape is a struct array of shapes which includes the shape data extracted
%    from the dbf file
%
% See also writeSHP, readExp, writeExp, plotShpP
%
% TO 090729


% Copyright 2009 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later


fprintf('\n\nplotshp --- %s\n',datestr(now));


%% Legal shape types
% shapetypes={
%     0, 'Null Shape'
%     1, 'Point'
%     3, 'PolyLine'
%     5, 'Polygon'
%     8, 'MultiPoint'
%     11,'PointZ'
%     13,'PolyLineZ'
%     15,'PolygonZ'
%     18,'MultiPointZ'
%     21,'PontM'
%     23,'PolyLineM'
%     25,'PolygonM'
%     28,'MultiPointM'
%     31,'MultiPatch'};

%% PartTypes defined for MultiPatch
% parttypes={
%     0,'TriangleStrip'
%     1,'TriangleFan'
%     2,'OuterRing'
%     3,'InnerRing'
%     4,'FirstRing'
%     5,'Ring'};

%% Skip the header because it is idential to that of the .shp file, where
%% we'll read it

for iShape=1:length(shape)
  
   switch shape(iShape).Type
       case 0  % Point, skip, only type was read
       case 1  % Point
         plot(shape(iShape).Points(:,1),shape(iShape).Points(:,2),'b.');
       case {3,5}  % 'PolyLine' 'Polygon'
           plotbox(shape(iShape).Box); hold on
           for iParts=1:length(shape(iShape).NumParts)
                first=shape(iShape).Parts(iParts)+1;
                if iParts<shape(iShape).NumParts
                    last=shape(iShape).Parts(iParts+1);
                else
                    last=shape(iShape).NumPoints;
                end
                range=first:last;
                plot(shape(iShape).Points(range,1),shape(iShape).Points(range,2));
            end
       case 8, % 'MultiPoint'
                plotbox=shape(iShape).Box;
                plot(shape(iShape).Points(:,1),shape(iShape).Points(:,2),'b.');
       case 11, % 'PointZ'  [x y z Measure]
          plot3(shape(iShape).Points(:,1),...
                shape(iShape).Points(:,2),...
                shape(iShape).Points(:,3),'b.'); hold on;
          text(shape(iShape).Points(:,1),...
                shape(iShape).Points(:,2),...
                shape(iShape).Points(:,3),sprintf('%s',...
                shape(iShape).Points(1,end)));
       case {13,15} % 'PolyLineZ' 'PolygonZ'
           plotbox(shape(iShape).Box,shape(iShape).ZAarray); hold on
           for iParts=1:length(shape(iShape).NumParts)
               first=shape(iShape).Parts(iParts)+1;
               if iParts<shape(iShape).NumParts
                   last=shape(iShape).Parts(iParts+1);
                else
                    last=shape(iShape).NumPoints;
                end
                range=first:last;
                plot3(shape(iShape).Points(range,1),...
                      shape(iShape).Points(range,2),...
                      shape(iShape).ZArray(range),'b');
                for iP=1:shape(iShape).NumPoints
                    text(shape(iShape).Points(iP,1),...
                         shape(iShape).Points(iP,2),...
                         shape(iShape).ZArray(iP),...
                         str2double(shape(iShape).MArray(iP)));
                end
            end
       case 18, % 'MultiPointZ'
           plotbox(shape(iShape).Box);
           for iP=1:shape(iShape).NumPoints
                plot3(shape(iShape).Points(range,1),...
                      shape(iShape).Points(range,2),...
                      shape(iShape).ZArray(range),'b.');
                text(shape(iShape).Points(iP,1),...
                    shape(iShape).Points(iP,2),...
                    shape(iShape).ZArray(iP),...
                    str2double(shape(iShape).MArray(iP)));
           end
       case 21, % 'PontM'          
           plot(shape(iShape).Points(:,1),shape(iShape).Points(:,2),'b.');
           text(shape(iShape).Points(:,1),shape(iShape).Points(:,2),str2double(shape(iShape).MArray(1)));
       case {23,25} % 'PolyLineM' and 'PolygonM'
           plotbox(shape(iShape).Box,shape(iShape).ZArray); hold on
           for iParts=1:length(shape(iShape).NumParts)
                first=shape(iShape).Parts(iParts)+1;
                if iParts<shape(iShape).NumParts
                    last=shape(iShape).Parts(iParts+1);
                else
                      
                end
                range=first:last;
                plot3(shape(iShape).Points(range,1),...
                      shape(iShape).Points(range,2),...
                      shape(iShape).ZArray(range  ),'b');
            end
           for iP=1:shape(iShape).NumPoints
                text(shape(iShape).Points(iP,1),...
                     shape(iShape).Points(iP,2),...
                     shape(iShape).ZArray(iP),...
                     str2double(shape(iShape).MArray(iP)));
           end
       case 28, % 'MultiPointM'
           plotbox(shape(iShape).Box,shape(iShape).ZArray); hold on
           for iP=1:shape(iShape).NumPoints
                plot3(shape(iShape).Points(range,1),...
                      shape(iShape).Points(range,2),...
                      shape(iShape).ZArray(range  ),'b.');
            end
           for iP=1:shape(iShape).NumPoints
                text(shape(iShape).Points(iP,1),...
                     shape(iShape).Points(iP,2),...
                     shape(iShape).ZArray(iP),...
                     str2double(shape(iShape).MArray(iP)));
           end
       case 31, % 'MultiPatch'};
           plotbox(shape(iShape).Box,shape(iShape).ZArray);
           for iP=1:shape(iShape).NumParts
               first=shape(iShape).Parts(iParts)+1;
               if iParts<shape(iShape).NumParts
                    last=shape(iShape).Parts(iParts+1);
               else
                    last=shape(iShape).NumPoints;
               end
               range=first:last;

               switch shape(iShape).PartType(iP)
                   case {0,1}  % triagle strip,  triangle fan
                        plot3(shape(iShape).NumPoints([1:3 1],1),...
                              shape(iShape).NumPoints([1:3 1],2),...
                              shape(iShape).ZArray([1:3 1]),'b');
                        for iT=4:shape(iShape).NumPoints
                            plot3(shape(iShape).NumPoints([iT-1 iT iT-2],1),...
                                  shape(iShape).NumPoints([iT-1 iT iT-2],2),...
                                  shape(iShape).ZArray(   [iT-1 iT iT-2]),'b');
                        end
                   case {2,3,4,5}  % Outer Ring, Inner Ring, First Ring, % Ring
                        plot3(shape(iShape).NumPoints(range,1),...
                              shape(iShape).NumPoints(range,2),...
                              shape(iShape).ZArray(range),'b');
               end
           end
       otherwise
           fprintf('Non implemented shape file %d skipped\n',shape.type)
   end
end

function plotbox(box,zrange)
if ~exist('zrange','var') || isempty(zrange) %2d
    plot(box(   [1 3 3 1 1]),...
         box(   [2 2 4 4 2]),...
         'r');  
else  % 3d    
plot3(box(   [1 3 3 1 1    1 3 3 1 1    1 3 3 1 1    1 3 3 1 1]),...
      box(   [2 2 4 4 2    2 2 4 4 2    2 2 2 2 2    4 4 4 4 4]),...
      zrange([1 1 1 1 1    2 2 2 2 2    1 1 2 2 1    1 1 2 2 1]),'r');
end

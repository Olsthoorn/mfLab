% mf_script2run showing use of getREGISIIXS(kmlpathfile,endDepth)
%
% REGIS is the national groundwater database of the Netherlands
% this function gets a cross section anywhere through the Netherlands
% directly from this database into your browser. The input file kmlfile
% is a pathfile you get from GooogleEarth. This pathfile is scanned, the
% coordinates are extracted and translated to the Dutch National Coordinate
% system and inserted into the URL, after which the cross section is
% requested from the REGIS databse. You can save the obtained page in your
% browser yourself directly from your browser. This yields a html file that
% will regenerate the cross section at any moment. Of course, you can also
% make a screendumpt if desired or a selection for your report.
%
% TO 100612

d=dir('*.kml');

%for i=1:length(d)
    i=1;
    
    fprintf('KML file is %s\n',d(i).name);
    getRegisIIXS(d(i).name);
    
%end


function [px py ix iy]=LatLonToPixels(N,E,zoom)
%  LatLonToPixels: Converts lat/lon to pixel coordinates in given zoom of the EPSG:4326 pyramid"

          fprintf('%15.10g',N);
          fprintf('%15.10g',E)

n=2^zoom;         fprintf('%15.10g',n);
w=1/n;            fprintf('%15.10g',w)
x=(E+180)/360;    fprintf('%15.10g',x);
y=(90-N)/180;     fprintf('%15.10g',y);
xn=x/w;           fprintf('%15.10g',xn);
yn=y/w;           fprintf('%15.10g',yn);
ix=fix(xn);       fprintf('%15.10g',ix);
iy=fix(yn);       fprintf('%15.10g',iy);
px =256*(xn-ix);  fprintf('%15.10g',px);
py =256*(yn-iy);  fprintf('%15.10g',py);
                  fprintf('\n');

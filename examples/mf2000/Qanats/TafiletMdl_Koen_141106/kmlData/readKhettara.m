name = 'khettara21thCentury';

[E,N]=kmlpath(name);
[x,y] = wgs2utm(N,E);
zDem = interp2(gr.xc,gr.yc,gr.Z(:,:,1),x,y);
zBot = interp2(gr.xc,gr.yc,gr.Z(:,:,end),x,y);

fprintf('%.9g\t%.9g\t%.9g\t%.9g\t%.5g\t%.5g\n',[N,E,x,y,zDem,zBot]');


%% An easy way to read objects and get them into Excel is as follows

%% example 1

p = khettara21thCentury.kml;
p.print;

% copy from command window directly to Excel


%% example 2

p = kmlPathsObj(fullfile('khettaras','khettaras.kml'));
p.print;

% copy from command window directly to Excel


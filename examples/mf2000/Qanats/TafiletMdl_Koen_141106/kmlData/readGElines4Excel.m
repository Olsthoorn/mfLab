% Hoe haal ik lijnen van GE binnen en hevel ze over naar Excel

% SEE ALSO:
% zie vooral kmlPathsObj dat object doet alles min of meer automatisch en
% berekent ook lengte en oppvervlak van de lijnen.
% which kmlPathsObj  (let op het meervoud, want werkt voor veel
% paden/lijnen tegelijk.
% help  kmlPathsObj  


toeHeels = {'Sayed' 16 1.5
            'Fougania' 15.5 1.5
            'Lakdima' 17.2 1.5
            'Bouia Lakadima' 17.1 1.5
            'Bouia Gadida' 15.9 1.5
            'Krayr Lakadima' 16.4 1.5
            'Mustaphia' 16.5 1.5
            'Aloria' 15 1.5
            'Krayr Gadida' 17.5 1.5
            'Grenia' 15 1.5
            'Auctania' 17.1 1.5
            };

if true
    p = kmlFolderObj(fullfile('khettaras','khettaras.kml'));

    p = p.getDem(gr);

    % voor khettaras
    depthToe = 17; depthHeel=1.5;
    p = p.setElevation(depthToe,depthHeel);

    p = p.setElevation(toeHeels);


    i=1; display(p(i).zDem-p(i).z) % zonder ; om te laten zien, diepte khettara (i) beneden dem


    p.print;
end

%% Read locations and dedpths objects

if true
    kmlFolderName = 'Super Khettara Fezna (straight).kml';
    % Hoe haal ik lijnen van GE binnen en hevel ze over naar Excel

    p = kmlFolderObj(kmlFolderName);

    p = p.getDem(gr);

    % voor khettaras
    depthToe = 17; depthHeel=1.5;
    p = p.setElevation(depthToe,depthHeel);

    p = p.setElevation(toeHeels);


    i=1; display(p(i).zDem-p(i).z); % zonder ; om te laten zien, diepte khettara (i) beneden dem


    p.print;
end

%% en copieer nu deze lijst van de matlab command window naar Excel
if false
    
    %% AREA OF A catchment (1)
    [E,N] = kmlpath('Wells.kml'); %#ok
    [xv,yv] = wgs2utm(N,E);

    % ===== What is het oppervlak binnen een contour?? ===========================
    % Dit kan je wiskundig precies brekenen, het gemakkelijkst met complex
    % getallen of met vector analyse.
    % Als je allemaal vectoren hebt vanaf een punt binnen een
    % polygong dan vormt elk stel vectoren een driehoek opgespannen door twee
    % vectoren vanaf dat centrale punt en een ribbe van het polygon. Als je
    % alle driehoeken optelt heb je het antwoord. Het bepalen van het oppervlak
    % van een driehoek gaat het gemakkelijk met het uitwendig product van twee
    % vectoren. Het anwoord daarvan is tweemaal dat opperlvak (plus of min
    % afhankelijk van de richting van de hoek ertussen.
    % extend x and y to a closed polygon

    %x = [0 1 1 0];
    %y = [0 0 1 1];

    xv=xv(:); yv=yv(:);
    if xv(end)~=xv(1) || yv(end)~=yv(1)
        xv = [xv; xv(1)];
        yv = [yv; yv(2)];
    end

    xp=mean(xv);
    yp=mean(yv);

    dx = xv-xp;
    dy = yv-yp;

    dA = NaN(3,numel(xv)-1);

    for i = size(dA,2):-1:1
        dA(:,i) =  cross([dx(i); dy(i); 0], [dx(i+1); dy(i+1); 0]);
    end
    A = sum(dA(3,:))/2;
    fprintf('A = %10g\n',A);


    %% AREA of a CATCHMENT 2
    % Methode twee. Kijk welke cellen van het grid binnen je polygon vallen en
    % sommer die en vergeminigvuldig dit met het cel oppervlak. Dan heb je het
    % ook met de nauwkeurighei van plus of min een halve ce.
    % Punt is dat een deel van het Catchment buiten het model gebied ligt.
    % We moeten daarom het hele DEM gebruiken voor de vergelijking.
    [DEMtot,E,N] = readASC(['..', filesep, 'DEM_SRTM', filesep, 'srtm_36_06.asc']);
    [xGr,yGr] = wgs2utm(N,E);
    xm = 0.5 * (xGr(1:end-1)+xGr(2:end));
    ym = 0.5 * (yGr(1:end-1)+yGr(2:end));
    dx =     diff(xGr);
    dy = abs(diff(yGr(:)));

    [Xm,Ym] = meshgrid(xm,ym);

    dA = dy*dx;

    I = inpolygon(Xm,Ym,x,y);
    nrOfCellsInI = sum(I(:));
    area = sum(sum(I .* dA));

    % Results
    fprintf('method1,   A = %10g\nmethod2,   A = %10g\n',A,area);
    fprintf('The difference is %.2g ha',abs(A-area)/1e4);


    plot(x,y,'y.-');
end

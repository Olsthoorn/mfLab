function tne = getMeteo(file)
    % Scan meteofile naar array A
    fid = fopen(file,'r');
    A   = fscanf(fid,'%d-%d-%d %f %f',[5,Inf])';
    fclose(fid); % inlezen klaar, sluit file

    % Splits de tijdreeks in [datum, neerslag en Makkinkverdamping]
    tne = [datenum(A(:,3),A(:,2),A(:,1)) A(:,[4 5])/1000]; % [t P N] % to mm/d

    % Papporteer de lengte van de ingelezen reeks
    fprintf('Length of time series = %d\n',length(tne(:,1)))
end

function [knmiData,knmiMeta] = readKNMI(knmiDataFile)
% imports KNMIdatafile. This function is obsolete
% see KNMIdata for help
%
% TO 141021

% knmiDataFile = 'KNMI_20131209.txt'; 

fid= fopen(knmiDataFile,'r');

fseek(fid,0,-1);

clear meta

knmiMeta{100,2} = '';

while true    
    s = fgetl(fid);
    n = regexp(s,'# STN','once');
    if ~isempty(n)
        break;
    end
end

k=1;
knmiMeta{k,1} = 'STN';

stn = fscanf(fid,'# %d:',1); %#ok
Lon = fscanf(fid,'%f',1);    %#ok
Lat = fscanf(fid,'%f',1);    %#ok
Alt = fscanf(fid,'%f',1);    %#ok

knmiMeta{k,2} = fscanf(fid,'%s',1);

fgetl(fid);
fgetl(fid);

while true
    
    s = fgetl(fid);
    n = regexp(s,'# [A-Z0-9]+ += ','once');
    if  ~isempty(n)
        k=k+1;
        knmiMeta{k,1} = s(regexp(s,'\w+ +=','once'):regexp(s,'\s+=','once')-1);
        knmiMeta{k,2} = s(regexp(s,'=','once')+2:end);                         
    else
        knmiMeta = knmiMeta(1:k,:);
        break;
    end
end

fgetl(fid);
fgetl(fid);

j=0;
knmiData = NaN(100000,size(knmiMeta,1));
fprintf('reading data from <<%s>> ...\n',knmiDataFile);
tic;
while true
    s = fgetl(fid);
    if s<0
        break;
    end
    j=j+1;
    if ~rem(j,50), fprintf('.'); end
    if ~rem(j,1000), fprintf('%d\n',j); end        
    knmiData(j,:) = cellfun(@strToNum,regexp(s,'[-0-9 ]+','match'));
end
knmiData = knmiData(1:j,:);

save KNMI knmiData knmiMeta

fprintf('%d\n',j);
fprintf('... done in %fs\ne',toc);


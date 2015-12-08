function [A names] = makeArray(xlsname,sheetnm,filesdir,columnHdr,VarNm,Ix,Iy)
%MAKEARRAY generates a 3D array of the size of the NHI model or part of it given xLim,yLim
%
% Example:
%    A = makeArray(gr,xlsname,sheetnam,filesdir,VarNm [,Ix [,Iy]])-- 
%
% TO 120401

warning('off') %#ok
[~,~,txt] = xlsread(xlsname,sheetnm,'','basic');
warning('on'); %#ok

hdr = txt(1    ,:);    % first line contains headers
txt = txt(2:end,:);    % values (char)

jVar         = strmatchi(columnHdr,hdr);
jZip         = strmatchi('zip'  ,hdr);
jFile        = strmatchi('file' ,hdr);
jDescription = strmatchi('descr',hdr);

I = strmatchi(VarNm,txt(:,jVar))';
[~,J] = sort(txt(I,jVar));
names{length(I),length(hdr)}='Dum';
for i=1:length(I)
    names(i,:) = {txt{I(J(i)),jVar} txt{I(J(i)),jFile} txt{I(J(i)),jZip} txt{I(J(i)),jDescription}};
end

%% Just try to find the next variable with the same basename
%  when it crashes we know the previos was the the last one
I=NaN(100,1); % 100 is just a dummy of sufficient size
for i=1:length(I)
    I(i)= strmatchi([VarNm sprintf('%d',i)],txt(:,jVar),'exact');   % col j contains file names
    if I(i)==0,
        Nlay=i-1;
        break;
    end
end
I=I(1:Nlay);

%% The number of layers is now know and also the locatin in the file overview
%  in the workbook
%  I is the row of the varialbe we seek, j its column, j+1 the column with
%  the corresponding file name (this is almost hard wiring but not too bad)

A=NaN(length(Iy),length(Ix),Nlay); % allocated

for iLay=1:Nlay
    A(:,:,iLay) = getNHIASC(fullfile(filesdir,txt{I(iLay),jFile}),Ix,Iy); % becomes a layer in 3D array
end
    
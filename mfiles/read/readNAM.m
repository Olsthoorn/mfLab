function nam=readNAM(fname,pth)
%READNAM reads the name file
%
% Example:
%    nam=readNAM(fname,pth)
%
% TO 070630 090713

% Copyright 2009-2013 Theo Olsthoorn, TU-Delft and Waternet, without any warranty
% under free software foundation GNU license version 3 or later

fprintf('# MATLAB readNAM %s\n',datestr(now));

fid=fopen([pth fname],'r');
nam=cell(100,3);
k=0;
while 1
    s=fgets(fid); if s==-1; break; end
    fprintf(s);
    if s(1)~='#' && ~isempty(sscanf(s,'%s'))  % skip comment and blank lines
        k=k+1;
        C=textscan(s,'%s %d %s',1);
        nam(k,1)=C{1};
        nam(k,2)=C(2);
        nam(k,3)=C{3};
    end
end
nam(k+1:end,:)=[];

%nam=textscan(fid,'%s'); nam=reshape(nam{1},3,length(nam{1})/3)';


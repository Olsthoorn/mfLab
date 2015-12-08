function [A]=readSCD(FName,pth,ptype,IxLim,IyLim)
%READSCD reads SCD (stress) files van NHI model (drains, rivers ghb and wells)
%
% Example:
%    A=readSCD(FName,pth,ptype,IxLim,IyLim)
%
% IxLim is 2-value vector holding lowest and highest col nr of AGV in NHI model
% IyLim is 2-value vector holdnig lowest and highers row nr of AGV in NHI model
%
% TO 120401

if pth(end)~='\', pth=[pth '\']; end
fid=fopen([pth FName],'r');

fgets(fid);
N=fscanf(fid,'%d',1);
fprintf('Reading SCD file %s, %d lines\n',FName,N);

switch ptype
    case 'd'  % drn_tot
        fmt='%10d%10d%10d%10f%15f';
        n=5;
    case 'w'  % wel_tot
        fmt='%10d%10d%10d%10f';
        n=4;
    case 'r'  % rivj
        fmt='%10d%10d%10d%10f%15f%10f%10f';
        n=7;
%     case 'r'  % rivw
%         fmt='%10d%10d%10d%10f%15f%10f%10f';
%         n=7;
%     case 'r'  % rivw_mz
%         fmt='%10d%10d%10d%10f%15f%10f%10f';
%         n=7;
     case 'g'
        fmt='%10d%10d%10d%10f%10f';
        n=5;
    otherwise
        error('Don''t know this type for SCD files ''%s''\n',ptype);
end


A=NaN(n,N);
for i=1:N
    v=fscanf(fid,fmt,[n,1]); fgets(fid);
    if v(3)>=IxLim(1) && v(3)<=IxLim(2) && v(2)>=IyLim(1) && v(2)<=IyLim(2)
       v(3)=v(3)-IxLim(1)+1;
       v(2)=v(2)-IyLim(1)+1;
       A(:,i)=v;
    end
end

A=A';
A=A(~isnan(A(:,1)),:);

fclose(fid);

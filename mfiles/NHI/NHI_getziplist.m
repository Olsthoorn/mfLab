function ziplist=NHI_getziplist(contents)
%NHI_GETZIPLIST  gets list of zip files contained in NHI on website
%   http://www.NHI.nu/downloads/
%   ziplist=NHI_zipfilelist(contents)
%
% INPUT:
%   contents is the name of the html files with the contents of the webiste
%     it can be obtained be downloading under right mouse button with the
%     browser open on the page showing the files.
%
%   Its feafult name is 
%      'Index of _downloads.html'
%   It uses this name of nargin=0;
%   The outcome is a celltarray with the names of the zip files on the page 
%
%   ziplist=NHI_getziplist;
%   ziplist=NHI_getziplist('Index_of_downloads.html');
%
% ToDo: thoroughly check (130429)
%
%   TO 110425

if nargin<1, contents='Index of _downloads.html'; end

fid=fopen(contents,'r');

ziplist=cell(300);

k=0;
while 1
    s=fgets(fid);
    if s==-1, break; end
    s=s(findstr('<a href=',s):end); s=s(findstr('>',s)+1:end); s=s(1:findstr(s,'</a>')-1);
    if findstr('.zip',s)>0,
        k=k+1;
        ziplist{k}=s;
    end
end
ziplist=ziplist(1:k);
fclose(fid);

fprintf('%d zip files on in file <<%s>>\n',length(ziplist),contents);

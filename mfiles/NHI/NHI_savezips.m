function NHI_savezips(ziplist,URL)
%NHI_SAVEZIPS saves the zip files of the NHI site
%
% Example:
%    matfile = NHI_savezips(ziplist[,URL])
%
%    ziplist obtained from NHI_getziplist
%    URL is URL of NHI site, default http://www.NHI.nu/downloads/
%    defaultmatfile = nhi.mat
%
%    See also: NHI_getziplist NHI_unzips
%
% ToDo: proper example
%
%    TO 110425

if nargin<2, URL='http://www.nhi.nu/downloads/'; end
if URL(end)~='/', URL=[URL '/']; end

for i=1:length(ziplist)
   [fname,status]=urlwrite([URL ziplist{i}],ziplist{i});
   if status~=1,
       fprintf('Failed to donwnload file <<%s>>\n',ziplist{i});
   else
       fprintf('Downloaded:  %s\n',fname);
   end
end



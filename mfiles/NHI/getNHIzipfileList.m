function ziplist=getNHIzipfileList(contents)
%GETNHIZIPFILELIST  gets list of zip files contained in NHI on website (www.NHI.nu)
%
%   NHI is the Dutch National Hydrologic Instrument available on the following website:
%   http://www.NHI.nu/downloads/
%
% Example:
%      ziplist=getNHIzipfileList(contents)
%
%   The contents is the name of the html file that contains the content of the
%   webiste. It can be obtained be downloading using the right mouse button
%   with the browser open on the page that shows the files.
%
%   Its default name is 
%   'Index of _downloads.html'
%   getNHIzipfileList uses this name if nargin==0;
%   The outcome is a cellt array with the names of the zip files on the page 
%
% Example:
%    ziplist=NHI_getziplist;
%    ziplist=NHI_getziplist('Index_of_downloads.html');
%
%  If it may not work, because the site has been changed. You can still use
%  this url and see what is in it. When on that site, just click on a file
%  to download it. This may be just as convenient.
%
%  Note that there may be more on the site than specified explicitly on the
%  page www.NHI.nu/Bibiliotheek. I found for instance the layer resistances
%  in it,(c_laag.zip), which was not referenced correctly in the mentioned
%  page. Also I found recharge.zip and startingheads.zip, which were not
%  references at all.
%
%   TO 110425 120428

URL = 'http://www.nhi.nu/downloads.html';
%if nargin<1, contents='Index of _downloads.html'; end

urlwrite(URL,'zipfilelist')
fid=fopen(contents,'r');

if fid<0, error('Can''t open URL <<%s>>.Check by hand.',contents); end

ziplist=cell(300);

k=0;
while 1
    s=fgets(fid);
    if s==-1, break; end
    s=s(strfind(s,'<a href='):end); s=s(strfind(s)+1:end,'>'); s=s(1:strfind(s,'</a>')-1);
    if strfind(s,'.zip')>0,
        k=k+1;
        ziplist{k}=s;
    end
end
ziplist=ziplist(1:k);
fclose(fid);

fprintf('%d zip files on in file <<%s>>\n',length(ziplist),contents);

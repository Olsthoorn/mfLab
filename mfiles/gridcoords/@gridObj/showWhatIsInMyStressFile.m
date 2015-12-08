function showWhatIsInMyStressFile(o,varargin)
%GRIDOBJ/PLOTWHATSINMYSTRESSFILE what is in the MODFLOW DRN intput file ?
% plots markers on figure corresponding to cells contained in stress file
% like ???.DNR, ???.GHB, ???.WEL

[ax   ,varargin] = getType(varargin,'axis',gca);
[FName,varargin] = getNext(varargin,'char','');


fid = fopen(FName);

if fid<0, error('%s: Can''t open file <<FName>>'); end

fgets(fid); fgets(fid); fgets(fid);

LRC = NaN(10000,3);    

k=0;

while k<size(LRC,1)
    s = fgets(fid);
    try
        LRC(k+1,:) = sscanf(s,'%d',[1,3]);
    catch %#ok
        break;
    end
    k=k+1;
end

if isempty(varargin), varargin{1}='ks'; end

if k<1,
    return;
    error('%s: No data in file <<%s>>',mfilename,FName); %#ok
end

LRC = LRC(1:k,:);
plot(ax,o.xm(LRC(:,3)),o.ym(LRC(:,2)),varargin{:});

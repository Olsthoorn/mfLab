function zipinfo=NHI_unzips(ziplist,URL,outdir,xLim,yLim)
%NHI_UNZIPS downloads, unzips and saves the NHI zip files URL
%
% Example:
%    matfile = NHI_savezips(ziplist[,URL,outdir,xLim,yLim])
%
%    ziplist obtained from NHI_getziplist
%    URL is URL of NHI site, default http://www.NHI.nu/downloads/
%    defaultmatfile = nhi.mat
%
%    See also: NHI_getziplist 
%
%    TO 110425

if nargin<5, xLim=[-Inf Inf]; yLim=[-Inf Inf]; end
if nargin<3
    outdir=[pwd '/temp/'];
    if ~exist(outdir,'dir'), mkdir(outdir); end
end
if nargin<2, URL='http://www.nhi.nu/downloads/'; end
if URL(end)~='/', URL=[URL '/']; end

%% Meta file for all data files
zipinfo=repmat(struct('zipname',''),[length(ziplist),1]);

%% Processing
for i=1:length(ziplist)
  zipinfo(i).zipname=ziplist{i};
  
  %% get the file from the web
  thisURL=[URL ziplist{i}];
  inzip  =unzip(thisURL,outdir);
  for j=1:length(inzip)
      [P,N,E]=fileparts(inzip{j});
      zipinfo(i).inzip{j}=[N E];
      zipinfo(i).ext{j}  =E;
  end
  
  %% move to outdir
  
  %% Get meta data
  try
      j=strmatchi('.met',zipinfo(i).ext);
      zipinfo(i).meta=NHI_readmeta([outdir zipinfo(i).inzip{j}]);

      %% get info about ascii files in outdir
  
      J=strmatchi('.asc',zipinfo(i).ext);
  
      % Get number of layers from file names
      NLAY=-Inf;
      for j=1:length(J)
          [Pth,Nm]=fileparts(zipinfo(i).inzip{J(j)});
          iLay=Nm(end)-'0'; if iLay>9 || iLay<1, iLay=NaN; end
          NLAY=max(NLAY,iLay);
      end
      zipinfo(i).NLAY=NLAY;

      %% get mate about ascii files from first one
      header=NHI_readASC(zipinfo(i).inzip{J(1)},outdir,xLim,yLim,'hdr');

      %NCOL   = header{strmatchi('ncols',       header(:,1)),2};
      %NROW   = header{strmatchi('nrows',       header(:,1)),2};
      %XLL    = header{strmatchi('xllcorner',   header(:,1)),2};
      %YLL    = header{strmatchi('yllcorner',   header(:,1)),2};
      %CELLSZ = header{strmatchi('cellsize',    header(:,1)),2};
      %NOVAL  = header{strmatchi('NODATA_value',header(:,1)),2};
      zipinfo(i).ASC_hdr=header;
  catch ME
      % just skip no ASC files in this directory
  end
    
end



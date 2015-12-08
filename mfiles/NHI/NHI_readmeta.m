function meta=NHI_readmeta(FNmeta)
%NHI_READMETA reads meta data belonging to NHI ascii data files (ESRI)
%
% Example:
%    meta=NH_READDATA(FNmeta)
%
%  meta is a struct having the field names specified in FNmeta
%
% Example meta data file NHI
% # Metadata NHI
% #
% # Algemeen:
% - Zip file          : onttrekkingen.zip
% - Naam              : Onttrekkingen --> PM (wordt verbeterd)
% - Type              : ASCII MODFLOW
% - Publicatie datum  : 21/01/2009
% - Versienr bestand  : v001
% - Versienr model    : v001
% - Producent         : NHI projectgroep
% 
% # Gedetailleerd:
% - Beschrijving      : Grondwater onttrekkingen per modellaag
% - Eenheid           : m^3/dag
% - Resolutie         : 250 m
% - Herkomst/Bron     : Provincaal Grondwaterregister aangevuld met regionale data
% - Legenda           : -
% - Procesbeschrijving: Toegekend aan modellagen op basis van beslisregels
% - Model catagorie   : Invoer
% - Model subcatagorie: Ondergrond
% - Deelrapport       : Onttrekkingen
%
%  TO 110425

fid=fopen(FNmeta,'r'); if fid<1, error('Can''t open file ''%s''',FNmeta); end

while 1
    s=fgets(fid);
    while s(1)=='#'
        s=fgets(fid);
    end  % skip blanks and comments
    if s(1)==-1,  break;    end
    
    i=findstr('-',s);
    j=findstr(':',s);
    
    if ~isempty(j)
        field=strtrim(s(i+1:j-1));
        field(findstr(field,' '))='_';
        field(findstr(field,'/'))='_';
        arg  =strtrim(s(j+2:end));
        arg(findstr( arg,' '))='_';
        meta.(field)=arg;
    end
end


fclose(fid);

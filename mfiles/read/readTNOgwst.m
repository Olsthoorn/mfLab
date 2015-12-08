function o = readTNOgwst(o,fname)
% TO 141027

% Titel:,,,,,,,,,,,
% Gebruikersnaam:,,,,,,,,,,,
% Periode aangevraagd:,01-01-1800,tot:,11-09-2014,,,,,,,,
% Gegevens beschikbaar:,27-01-1961,tot:,31-12-2010,,,,,,,,
% Datum: ,11-09-2014,,,,,,,,,,
% Referentie:,NAP,,,,,,,,,,
% 
% NAP:,Normaal Amsterdams Peil,,,,,,,,,,
% MV:,Maaiveld,,,,,,,,,,
% MP:,Meetpunt,,,,,,,,,,
% 
% Locatie,Filternummer,Externe aanduiding,X-coordinaat,Y-coordinaat,Maaiveld (cm t.o.v. NAP),Datum maaiveld gemeten,Startdatum,Einddatum,Meetpunt (cm t.o.v. NAP),Meetpunt (cm t.o.v. MV),Bovenkant filter (cm t.o.v. NAP),Onderkant filter (cm t.o.v. NAP)
% B06D0110,001,06DP0110,190573,578183,70,13-07-1971,25-07-1996,14-01-2003,29,-41,-1293,-1393
% B06D0110,001,06DP0110,190573,578183,55,14-01-2003,14-01-2003,31-12-2010,49,-6,-1273,-1373
% 
% 
% Locatie,Filternummer,Peildatum,Stand (cm t.o.v. MP),Stand (cm t.o.v. MV),Stand (cm t.o.v. NAP),Bijzonderheid,Opmerking,,,
% B06D0110,001,27-01-1961,118,,,,,,,,
% B06D0110,001,27-02-1961,114,,,,,,,,
% B06D0110,001,27-03-1961,125,,,,,,,,

fp =  fopen(fname,'r');

% Locatie
% Filternummber
% Externe_aanduiding
% xcoord      =
% ycoord      =
% mvNAP       =
% mvdatum     = datenum(C{2}{ 9},C{2}{ 8},C{2}{ 7});
% startdatum  = datenum(C{2}{12},C{2}{11},C{2}{10});
% einddatem   = datenum(C{2}{15},C{2}{14},C{2}{13});
% mpntNAP     =
% mpntMV      =
% bkfNAP      =
% okfNAP      =

o.comment    = fname;
o.type       = 'head';

  date       = lookfor(fp,'Datum');
o.date       = ['TNO file date= ' date{end}];

Hdr          = lookfor(fp,'Locatie'); %#ok

C   = textscan(fp,'%s %d %s %f %f %f %f-%f-%f %f-%f-%f %f-%f-%f %f %f %f %f',2,'delimiter',',','emptyValue',NaN);

o.tnocode    = C{1}{end};
o.filtnr     = C{2}(end);

o.name       = sprintf('%s_%d',o.tnocode,o.filternr);

o.xcoord     = C{4}(end);
o.ycoord     = C{5}(end);

o.surflev    = C{16}(end)/100;
o.upfiltlev  = C{18}(end)/100;
o.lowfiltlev = C{19}(end)/100;

%% Non-Menyanthes info --> UserData

o.UserData.extName      = C{3}{end};
o.UserData.mvDatum      = datenum(C{ 9},C{ 8},C{ 7});
o.UserData.startDate    = datenum(C{12},C{11},C{10});
o.UserData.endDate      = datenum(C{15},C{14},C{13});
o.UserData.mpntNAP      = C{16}/100;
o.UserData.mpntMV       = C{17}/100;
o.UserData.bkfNAP       = C{18}/100;
o.UserData.okfNAP       = C{19}/100;

%% Locatie,Filternummer,Peildatum,Stand (cm t.o.v. MP),Stand (cm t.o.v. MV),Stand (cm t.o.v. NAP),Bijzonderheid,Opmerking,,,

  Hdr = lookfor(fp,'Locatie');
o.UserData.DataHdr = Hdr(2:end);

% Data
C = textscan(fp,'%s %f %f-%f-%f %f %f %f %s %s %f %f %f %f %f',Inf,'delimiter',',','emptyValue',NaN);

o.values = [datenum(C{5},C{4},C{3}),  [C{6}, C{7}, C{8}]/100 ];

o.UserData.remark      = {o.values(:,1) C{9} C{10}};

%fprintf('hello world\n');

end

function tokens = lookfor(fp,word)

for i=1:100
    s = fgetl(fp);
    n = length(word);
    if length(s)>=n && strcmpi(s(1:n),word)
        [tokens] = regexp(s, '[A-Za-z0-9\_\+\-\.\(\) ]*', 'match');         
        return;        
    end
end

end

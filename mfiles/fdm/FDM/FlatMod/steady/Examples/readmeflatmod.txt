% Readmeflatmod
% PAO cursus rekenen zonder kant- en klare modellen
% toelichting op de directory flatmod

Deze directory bevat de voorbeelden van numeriek rekenen met vlakke eindige differentiemodellen
TO 000820


flatmeshctrd.m		Vlak mesh-centered stationair eindig differentiemodel	(Als eindige elemente, heads in knopen)
flatblockctrd.m	Vlak blok-centered stationair eindig differentiemodel (soort Modflow, cellen)


idkBethune.xls		% doorsnede Bethune Wassen (1991) geeft de zonenummers aan van elke cel in doorsnede
idkBethune.txt		% zelfde als idBethune.xls, in txt formaat om ingelezen te worden bij voorbeeld


nmzex1.m
nmzex1commented.m = als nmzex1 met commentaar en pauzes, zodat berekening in stappen wordt uitgevoerd

nmzex1 berekent met flatmeshctrd een drielagenmodel uit. De scheidende lagen worden gemaakt door cellen
een kleine doorlatendheid te geven. In de tekening worden deze paars gekleurd. Tevens wordt de n-lagen
analytische oplossing berekend en als lijnen weergegeven ter vergelijking met het numerieke resultaat (plusjes).
In de onderset grafiek van figuur 1 wordt de analytisch berkende flux [m2/d] weergegeven die is berekend
in de watervoerende pakketten.
Figuur twee geeft de stijghoogtelijnen en de stroomlijnen weer. De stroomlijnen zijn verkregen door het zogenoemde
toegevoegde probleem uit te rekeken, dit levert de stroomfunctie. Hiertoe wordt met hetzelfde model nogmaals
gedraaid, nadat in elke cel de kx is vervangen door 1/ky en ky door 1/kx en de randvoorwaarden zijn aangepast. De
randvoorwaarde van het toegevoegde probleem is de totale stroming over de rand vanaf een bepaal punt op de rand.
Deze stroming wordt gehaald uit de oplossing van het stijghoogteprobleem. Bij mesh-centered modellen is het
altijd lastig om precies uit te maken hoeveel water tussen twee knoopppunten op de modelrand stroomt, omdat
alleen stijghoogten en volumestromen in de knooppunten worden berekend. Bij een block-centered model is het
omgekeerde het geval. De volumestroom is daar juist tussen de knooppunten (dus over de celwanden) bekend.
De derge grafiek berekent de cumulatiev verblijftijdsverdeling. Dit gebeurt langs de zojuist berekende stroomlijnen.
We weten immers de hoeveelheid water tussen elk paar stroomlijnen. We weten ook dat het water precies langs de
contourlijnen van de stroomfunctie stroomt. We kunnen nu langs deze stroomlijnen de snelheid in elk punt berekenen
en dus ook de verblijftijd. Door de verblijftijden van alle stroomlijnen acheraf bij elkaar te voegen, wordt
de cumulatieve verblijftijdsverdeling verkregen.

nmzex2.m     (te groot voor de student edition)     
nmzex2EDU.m  (Student edition, zie verderop, gebruikt kleinere dataset (zie EDU in naam)

Hier wordt het model flatmeshctrd toegepast op een veel grotere doorsnede, namelijk een die loopt van de Loos-
drechtse plassen door de Bethunepolder en dan naar de Utrechtse Heuvelrug (Wassen, 1991). Ook hier worden
stijghoogte- en stroomlijnen berekend en tenslotte de verblijftijdsverdeling.
De verdeling van de doorlatendheid over de cellen van de doorsnede wordt ingelezen uit een tekstbestand die
in een spreadsheet is aangemaakt. Deze matrix wordt ook gebruikt om de verschillende doorlatendheiden in de
doorsnede te tekenen.

nmzex2EDU.m

De dataset die door nmzex2 wordt gebruikt is te groot voor de studentenversie van MATLAB!
Speciaal voor de Student Edition van Matlab is hier een kleinere versie (dataset) die
wel in de Student Edition kan worden gedraaid. Ik heb de dataset grofweg verkleind, zonder heel precies
te letten op de consequenties die dit heeft voor de randvoorwaarden. Dus kleine afwijkijkingen kunnen
hierdoor tussen het grotere en het kleinere model optreden.


MeshvsBlock.m		% voorbeeld waarin resultaten mesh- en block-centered model worden vergeleken.

Hier vergelijken we het mesh-gecentreerde model met het block-gecentreerde (MODFLOW-achtige) model.
Na elk van beide berekeningen pauzeert de mfile en moet op ENTER worden gedrukt om door te gaan.
Het mesh-gecentreerde model berekent stijghoogten in knooppunten. Dit zijn de hoeken van de cellen of
elementen. Dit is boven- en onderin elke cellenlaag, waardoor goede contourlijnen worden verkregen. Het
blok-gecentreerde model berekent stijghoogten alleen in het midden van de cellen. Contouren van een doorsnede
zoals hier, waarin goed- en slechtdoorlatende lagen elkaar afwisselen komen dan niet overeen met de werkelijkheid,
omdat het conoutouring programma alleen kijkt naar het midden van de cellen. Dat de resultaten van beide
modellen wel met elkaar in overeenstemming zijn blijkt uit het feit dat de contouren van beide modellen elkaar
in het midden van de cellen precies snijden. Ook de vergelijking met de analytische meer-lagen oplossing in
figuur 2 laat zien dat de resultaten gelijk zijn.
Deze verkeerde contouren krijg je dus altijd als je dwarsdoornsneden weergeeft van een blok-gecentreerde model,
dat een opeenstapeling is van WVP's en SDP's. Dit is bijvoorbeeld standaard het geval bij Visual Modflow. Je
kan dit oplossen door elk WVP te splitsen in 3 afzonderlijke cellenlagen, waarbij de laag boven- en onderin
het WVP zeer dun zijn. In MATLAB kun je dit ook manipuleren door de stijghoogte die in het midden van WVP's
worden berekend toe te kennen aan de boven- en onderzijde van de WVP en dan pas te contouren.

sysmat.m          systeemmatrix nodig voor analytische n-lagen berkening
getcontours.m		ontrafeld contouren om er de verblijftijd mee uit te rekenen
pauzeer.m			pauzeren tot gebruiker op ENTER drukt

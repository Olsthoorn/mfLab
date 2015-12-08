Beide modellen zijn stationaire, 2D blokgecentreerde eindige differentiemodellen die behalve de stijghoogte (in de cellen) ook de stroomfunctie berekenen (in de knooppunten). Dit gebeurt met de flows over de celwanden die bekend zijn zodra de stijghoogten zijn berekend. Het radiale en vlakkemodel zijn precies gelijk behoudens de wijze waarop de matrixcoefficienten zijn berekend. Het radiale model houdt uiteraard rekening met de toenemende omvang van de ringen. Het vlakke model heeft dit niet.

De voorbeelden zijn bijgevoegd. De stroomfuntie kan gebruikt of genegeerd worden. Er worden altijd twee stroomfuncties berekend, een met verticale branchcuts en een met horizontale. De te gebruiken versie hangt af van de wijze waarop watet van interne bronnen wordt afgevoerd (naar horzonale bovendrand of verticale zijrand). Zonder interne bronnen zijn beide stroomfunctiematrices aan elkaar gelijk.

[Phi,Q,PsiY,PsiX]=flatblockctrd(x,r,Kx,Kr,FH,Q);
[Phi,Q,PhiY,PhiX]= radblockctrd(x,z,Kx,Kz,FH,Q);

TO 000904
Read("CLP_core.g");;  Read("pairs_mid.g");;  P1 := PAIRS;;
Read("pairs_bigm.g");;  PAIRS := Concatenation( P1, PAIRS );;
pr := First( PAIRS, p -> p.l = 108 and p.i = 16 and p.m = 54 and p.lift = 18 );;
P := LoadPair( pr );;  PrepareActions( P );;  Print( P.desc, "\n" );
for t in [ 1 .. 300 ] do
A := [ [ SupportWithOne(P.nA,Random([1,2])), SupportWithOne(P.nA,Random([1,2])) ], [ SupportWithOne(P.nA,Random([1,2])), RandomEntry(P.nA,1,2) ] ];;
B := [ [ SupportWithOne( P.nL, 3 ) ] ];;
if ForAny( Flat(A), x -> not P.AinN[x] ) then
c := CLP_Build( P, A, B, [2,2], [1,1] );;
if not HasCancellation( c, A, B, [2,2], [1,1], P.m ) then
K := c.N - RankMat(c.HX) - RankMat(c.HZ);;
if K >= 8 then
dZ := DistRandCSS( c.HX, c.HZ, 300, 0 : field := GF(2) );;
dX := DistRandCSS( c.HZ, c.HX, 300, 0 : field := GF(2) );;
Print( "K=", K, "  rankRb=", RankMat( c.RB[1][1] ), "/", P.m, "  dZ<=", dZ, "  dX<=", dX, "\n" );
fi; fi; fi;
od;
QUIT;
